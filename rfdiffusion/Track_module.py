import torch.utils.checkpoint as checkpoint
from rfdiffusion.util_module import *
from rfdiffusion.Attention_module import *
from rfdiffusion.SE3_network import SE3TransformerWrapper

# Components for three-track blocks
# 1. MSA -> MSA update (biased attention. bias from pair & structure)
# 2. Pair -> Pair update (biased attention. bias from structure)
# 3. MSA -> Pair update (extract coevolution signal)
# 4. Str -> Str update (node from MSA, edge from Pair)

# Update MSA with biased self-attention. bias from Pair & Str
class MSAPairStr2MSA(nn.Module):
    def __init__(self, d_msa=256, d_pair=128, n_head=8, d_state=16,
                 d_hidden=32, p_drop=0.15, use_global_attn=False):
        super(MSAPairStr2MSA, self).__init__()
        self.norm_pair = nn.LayerNorm(d_pair)
        self.proj_pair = nn.Linear(d_pair+36, d_pair)
        self.norm_state = nn.LayerNorm(d_state)
        self.proj_state = nn.Linear(d_state, d_msa)
        self.drop_row = Dropout(broadcast_dim=1, p_drop=p_drop)
        self.row_attn = MSARowAttentionWithBias(d_msa=d_msa, d_pair=d_pair,
                                                n_head=n_head, d_hidden=d_hidden) 
        if use_global_attn:
            self.col_attn = MSAColGlobalAttention(d_msa=d_msa, n_head=n_head, d_hidden=d_hidden) 
        else:
            self.col_attn = MSAColAttention(d_msa=d_msa, n_head=n_head, d_hidden=d_hidden) 
        self.ff = FeedForwardLayer(d_msa, 4, p_drop=p_drop)
        
        # Do proper initialization
        self.reset_parameter()

    def reset_parameter(self):
        # initialize weights to normal distrib
        self.proj_pair = init_lecun_normal(self.proj_pair)
        self.proj_state = init_lecun_normal(self.proj_state)

        # initialize bias to zeros
        nn.init.zeros_(self.proj_pair.bias)
        nn.init.zeros_(self.proj_state.bias)

    def forward(self, msa, pair, rbf_feat, state):
        '''
        Inputs:
            - msa: MSA feature (B, N, L, d_msa)
            - pair: Pair feature (B, L, L, d_pair)
            - rbf_feat: Ca-Ca distance feature calculated from xyz coordinates (B, L, L, 36)
            - xyz: xyz coordinates (B, L, n_atom, 3)
            - state: updated node features after SE(3)-Transformer layer (B, L, d_state)
        Output:
            - msa: Updated MSA feature (B, N, L, d_msa)
        '''
        B, N, L = msa.shape[:3]

        # prepare input bias feature by combining pair & coordinate info
        pair = self.norm_pair(pair)
        pair = torch.cat((pair, rbf_feat), dim=-1)
        pair = self.proj_pair(pair) # (B, L, L, d_pair)
        #
        # update query sequence feature (first sequence in the MSA) with feedbacks (state) from SE3
        state = self.norm_state(state)
        state = self.proj_state(state).reshape(B, 1, L, -1)
        msa = msa.index_add(1, torch.tensor([0,], device=state.device), state)
        #
        # Apply row/column attention to msa & transform 
        msa = msa + self.drop_row(self.row_attn(msa, pair))
        msa = msa + self.col_attn(msa)
        msa = msa + self.ff(msa)

        return msa

class PairStr2Pair(nn.Module):
    def __init__(self, d_pair=128, n_head=4, d_hidden=32, d_rbf=36, p_drop=0.15):
        super(PairStr2Pair, self).__init__()
        
        self.emb_rbf = nn.Linear(d_rbf, d_hidden)
        self.proj_rbf = nn.Linear(d_hidden, d_pair)

        self.drop_row = Dropout(broadcast_dim=1, p_drop=p_drop)
        self.drop_col = Dropout(broadcast_dim=2, p_drop=p_drop)

        self.row_attn = BiasedAxialAttention(d_pair, d_pair, n_head, d_hidden, p_drop=p_drop, is_row=True)
        self.col_attn = BiasedAxialAttention(d_pair, d_pair, n_head, d_hidden, p_drop=p_drop, is_row=False)

        self.ff = FeedForwardLayer(d_pair, 2)
        
        self.reset_parameter()
    
    def reset_parameter(self):
        nn.init.kaiming_normal_(self.emb_rbf.weight, nonlinearity='relu')
        nn.init.zeros_(self.emb_rbf.bias)
        
        self.proj_rbf = init_lecun_normal(self.proj_rbf)
        nn.init.zeros_(self.proj_rbf.bias)

    def forward(self, pair, rbf_feat):
        B, L = pair.shape[:2]

        rbf_feat = self.proj_rbf(F.relu_(self.emb_rbf(rbf_feat)))

        pair = pair + self.drop_row(self.row_attn(pair, rbf_feat))
        pair = pair + self.drop_col(self.col_attn(pair, rbf_feat))
        pair = pair + self.ff(pair)
        return pair

class MSA2Pair(nn.Module):
    def __init__(self, d_msa=256, d_pair=128, d_hidden=32, p_drop=0.15):
        super(MSA2Pair, self).__init__()
        self.norm = nn.LayerNorm(d_msa)
        self.proj_left = nn.Linear(d_msa, d_hidden)
        self.proj_right = nn.Linear(d_msa, d_hidden)
        self.proj_out = nn.Linear(d_hidden*d_hidden, d_pair)
        
        self.reset_parameter()

    def reset_parameter(self):
        # normal initialization
        self.proj_left = init_lecun_normal(self.proj_left)
        self.proj_right = init_lecun_normal(self.proj_right)
        nn.init.zeros_(self.proj_left.bias)
        nn.init.zeros_(self.proj_right.bias)

        # zero initialize output
        nn.init.zeros_(self.proj_out.weight)
        nn.init.zeros_(self.proj_out.bias)

    def forward(self, msa, pair):
        B, N, L = msa.shape[:3]
        msa = self.norm(msa)
        left = self.proj_left(msa)
        right = self.proj_right(msa)
        right = right / float(N)
        out = einsum('bsli,bsmj->blmij', left, right).reshape(B, L, L, -1)
        out = self.proj_out(out)
       
        pair = pair + out
        
        return pair

class SCPred(nn.Module):
    def __init__(self, d_msa=256, d_state=32, d_hidden=128, p_drop=0.15):
        super(SCPred, self).__init__()
        self.norm_s0 = nn.LayerNorm(d_msa)
        self.norm_si = nn.LayerNorm(d_state)
        self.linear_s0 = nn.Linear(d_msa, d_hidden)
        self.linear_si = nn.Linear(d_state, d_hidden)

        # ResNet layers
        self.linear_1 = nn.Linear(d_hidden, d_hidden)
        self.linear_2 = nn.Linear(d_hidden, d_hidden)
        self.linear_3 = nn.Linear(d_hidden, d_hidden)
        self.linear_4 = nn.Linear(d_hidden, d_hidden)

        # Final outputs
        self.linear_out = nn.Linear(d_hidden, 20)

        self.reset_parameter()

    def reset_parameter(self):
        # normal initialization
        self.linear_s0 = init_lecun_normal(self.linear_s0)
        self.linear_si = init_lecun_normal(self.linear_si)
        self.linear_out = init_lecun_normal(self.linear_out)
        nn.init.zeros_(self.linear_s0.bias)
        nn.init.zeros_(self.linear_si.bias)
        nn.init.zeros_(self.linear_out.bias)
        
        # right before relu activation: He initializer (kaiming normal)
        nn.init.kaiming_normal_(self.linear_1.weight, nonlinearity='relu')
        nn.init.zeros_(self.linear_1.bias)
        nn.init.kaiming_normal_(self.linear_3.weight, nonlinearity='relu')
        nn.init.zeros_(self.linear_3.bias)

        # right before residual connection: zero initialize
        nn.init.zeros_(self.linear_2.weight)
        nn.init.zeros_(self.linear_2.bias)
        nn.init.zeros_(self.linear_4.weight)
        nn.init.zeros_(self.linear_4.bias)
    
    def forward(self, seq, state):
        '''
        Predict side-chain torsion angles along with backbone torsions
        Inputs:
            - seq: hidden embeddings corresponding to query sequence (B, L, d_msa)
            - state: state feature (output l0 feature) from previous SE3 layer (B, L, d_state)
        Outputs:
            - si: predicted torsion angles (phi, psi, omega, chi1~4 with cos/sin, Cb bend, Cb twist, CG) (B, L, 10, 2)
        '''
        B, L = seq.shape[:2]
        seq = self.norm_s0(seq)
        state = self.norm_si(state)
        si = self.linear_s0(seq) + self.linear_si(state)

        si = si + self.linear_2(F.relu_(self.linear_1(F.relu_(si))))
        si = si + self.linear_4(F.relu_(self.linear_3(F.relu_(si))))

        si = self.linear_out(F.relu_(si))
        return si.view(B, L, 10, 2)


class Str2Str(nn.Module):
    def __init__(self, d_msa=256, d_pair=128, d_state=16, 
            SE3_param={'l0_in_features':32, 'l0_out_features':16, 'num_edge_features':32}, p_drop=0.1):
        super(Str2Str, self).__init__()
        
        # initial node & pair feature process
        self.norm_msa = nn.LayerNorm(d_msa)
        self.norm_pair = nn.LayerNorm(d_pair)
        self.norm_state = nn.LayerNorm(d_state)
    
        self.embed_x = nn.Linear(d_msa+d_state, SE3_param['l0_in_features'])
        self.embed_e1 = nn.Linear(d_pair, SE3_param['num_edge_features'])
        self.embed_e2 = nn.Linear(SE3_param['num_edge_features']+36+1, SE3_param['num_edge_features'])
        
        self.norm_node = nn.LayerNorm(SE3_param['l0_in_features'])
        self.norm_edge1 = nn.LayerNorm(SE3_param['num_edge_features'])
        self.norm_edge2 = nn.LayerNorm(SE3_param['num_edge_features'])
        
        self.se3 = SE3TransformerWrapper(**SE3_param)
        self.sc_predictor = SCPred(d_msa=d_msa, d_state=SE3_param['l0_out_features'],
                                   p_drop=p_drop)
        
        self.reset_parameter()

    def reset_parameter(self):
        # initialize weights to normal distribution
        self.embed_x = init_lecun_normal(self.embed_x)
        self.embed_e1 = init_lecun_normal(self.embed_e1)
        self.embed_e2 = init_lecun_normal(self.embed_e2)

        # initialize bias to zeros
        nn.init.zeros_(self.embed_x.bias)
        nn.init.zeros_(self.embed_e1.bias)
        nn.init.zeros_(self.embed_e2.bias)
    
    @torch.cuda.amp.autocast(enabled=False)
    def forward(self, msa, pair, R_in, T_in, xyz, state, idx, motif_mask, top_k=64, eps=1e-5):
        B, N, L = msa.shape[:3]

        if motif_mask is None:
            motif_mask = torch.zeros(L).bool()
        
        # process msa & pair features
        node = self.norm_msa(msa[:,0])
        pair = self.norm_pair(pair)
        state = self.norm_state(state)
       
        node = torch.cat((node, state), dim=-1)
        node = self.norm_node(self.embed_x(node))
        pair = self.norm_edge1(self.embed_e1(pair))
        
        neighbor = get_seqsep(idx)
        rbf_feat = rbf(torch.cdist(xyz[:,:,1], xyz[:,:,1]))
        pair = torch.cat((pair, rbf_feat, neighbor), dim=-1)
        pair = self.norm_edge2(self.embed_e2(pair))
        
        # define graph
        if top_k != 0:
            G, edge_feats = make_topk_graph(xyz[:,:,1,:], pair, idx, top_k=top_k)
        else:
            G, edge_feats = make_full_graph(xyz[:,:,1,:], pair, idx, top_k=top_k)
        l1_feats = xyz - xyz[:,:,1,:].unsqueeze(2)
        l1_feats = l1_feats.reshape(B*L, -1, 3)
        
        # apply SE(3) Transformer & update coordinates
        shift = self.se3(G, node.reshape(B*L, -1, 1), l1_feats, edge_feats)

        state = shift['0'].reshape(B, L, -1) # (B, L, C)
        
        offset = shift['1'].reshape(B, L, 2, 3)
        offset[:,motif_mask,...] = 0            # NOTE: motif mask is all zeros if not freeezing the motif 

        delTi = offset[:,:,0,:] / 10.0 # translation
        R = offset[:,:,1,:] / 100.0 # rotation
        
        Qnorm = torch.sqrt( 1 + torch.sum(R*R, dim=-1) )
        qA, qB, qC, qD = 1/Qnorm, R[:,:,0]/Qnorm, R[:,:,1]/Qnorm, R[:,:,2]/Qnorm

        delRi = torch.zeros((B,L,3,3), device=xyz.device)
        delRi[:,:,0,0] = qA*qA+qB*qB-qC*qC-qD*qD
        delRi[:,:,0,1] = 2*qB*qC - 2*qA*qD
        delRi[:,:,0,2] = 2*qB*qD + 2*qA*qC
        delRi[:,:,1,0] = 2*qB*qC + 2*qA*qD
        delRi[:,:,1,1] = qA*qA-qB*qB+qC*qC-qD*qD
        delRi[:,:,1,2] = 2*qC*qD - 2*qA*qB
        delRi[:,:,2,0] = 2*qB*qD - 2*qA*qC
        delRi[:,:,2,1] = 2*qC*qD + 2*qA*qB
        delRi[:,:,2,2] = qA*qA-qB*qB-qC*qC+qD*qD

        Ri = einsum('bnij,bnjk->bnik', delRi, R_in)
        Ti = delTi + T_in #einsum('bnij,bnj->bni', delRi, T_in) + delTi
            
        alpha = self.sc_predictor(msa[:,0], state)
        return Ri, Ti, state, alpha

class IterBlock(nn.Module):
    def __init__(self, d_msa=256, d_pair=128,
                 n_head_msa=8, n_head_pair=4,
                 use_global_attn=False,
                 d_hidden=32, d_hidden_msa=None, p_drop=0.15,
                 SE3_param={'l0_in_features':32, 'l0_out_features':16, 'num_edge_features':32}):
        super(IterBlock, self).__init__()
        if d_hidden_msa == None:
            d_hidden_msa = d_hidden

        self.msa2msa = MSAPairStr2MSA(d_msa=d_msa, d_pair=d_pair,
                                      n_head=n_head_msa,
                                      d_state=SE3_param['l0_out_features'],
                                      use_global_attn=use_global_attn,
                                      d_hidden=d_hidden_msa, p_drop=p_drop)
        self.msa2pair = MSA2Pair(d_msa=d_msa, d_pair=d_pair,
                                 d_hidden=d_hidden//2, p_drop=p_drop)
                                 #d_hidden=d_hidden, p_drop=p_drop)
        self.pair2pair = PairStr2Pair(d_pair=d_pair, n_head=n_head_pair, 
                                      d_hidden=d_hidden, p_drop=p_drop)
        self.str2str = Str2Str(d_msa=d_msa, d_pair=d_pair,
                               d_state=SE3_param['l0_out_features'],
                               SE3_param=SE3_param,
                               p_drop=p_drop)

    def forward(self, msa, pair, R_in, T_in, xyz, state, idx, motif_mask, use_checkpoint=False):
        rbf_feat = rbf(torch.cdist(xyz[:,:,1,:], xyz[:,:,1,:]))
        if use_checkpoint:
            msa = checkpoint.checkpoint(create_custom_forward(self.msa2msa), msa, pair, rbf_feat, state)
            pair = checkpoint.checkpoint(create_custom_forward(self.msa2pair), msa, pair)
            pair = checkpoint.checkpoint(create_custom_forward(self.pair2pair), pair, rbf_feat)
            R, T, state, alpha = checkpoint.checkpoint(create_custom_forward(self.str2str, top_k=0), msa, pair, R_in, T_in, xyz, state, idx, motif_mask)
        else:
            msa = self.msa2msa(msa, pair, rbf_feat, state)
            pair = self.msa2pair(msa, pair)
            pair = self.pair2pair(pair, rbf_feat)
            R, T, state, alpha = self.str2str(msa, pair, R_in, T_in, xyz, state, idx, motif_mask=motif_mask, top_k=0) 
        
        return msa, pair, R, T, state, alpha

class IterativeSimulator(nn.Module):
    def __init__(self, n_extra_block=4, n_main_block=12, n_ref_block=4,
                 d_msa=256, d_msa_full=64, d_pair=128, d_hidden=32,
                 n_head_msa=8, n_head_pair=4,
                 SE3_param_full={'l0_in_features':32, 'l0_out_features':16, 'num_edge_features':32},
                 SE3_param_topk={'l0_in_features':32, 'l0_out_features':16, 'num_edge_features':32},
                 p_drop=0.15):
        super(IterativeSimulator, self).__init__()
        self.n_extra_block = n_extra_block
        self.n_main_block = n_main_block
        self.n_ref_block = n_ref_block
        
        self.proj_state = nn.Linear(SE3_param_topk['l0_out_features'], SE3_param_full['l0_out_features'])
        # Update with extra sequences
        if n_extra_block > 0:
            self.extra_block = nn.ModuleList([IterBlock(d_msa=d_msa_full, d_pair=d_pair,
                                                        n_head_msa=n_head_msa,
                                                        n_head_pair=n_head_pair,
                                                        d_hidden_msa=8,
                                                        d_hidden=d_hidden,
                                                        p_drop=p_drop,
                                                        use_global_attn=True,
                                                        SE3_param=SE3_param_full)
                                                        for i in range(n_extra_block)])

        # Update with seed sequences
        if n_main_block > 0:
            self.main_block = nn.ModuleList([IterBlock(d_msa=d_msa, d_pair=d_pair,
                                                       n_head_msa=n_head_msa,
                                                       n_head_pair=n_head_pair,
                                                       d_hidden=d_hidden,
                                                       p_drop=p_drop,
                                                       use_global_attn=False,
                                                       SE3_param=SE3_param_full)
                                                       for i in range(n_main_block)])

        self.proj_state2 = nn.Linear(SE3_param_full['l0_out_features'], SE3_param_topk['l0_out_features'])
        # Final SE(3) refinement
        if n_ref_block > 0:
            self.str_refiner = Str2Str(d_msa=d_msa, d_pair=d_pair,
                                       d_state=SE3_param_topk['l0_out_features'],
                                       SE3_param=SE3_param_topk,
                                       p_drop=p_drop)
    
        self.reset_parameter()
    def reset_parameter(self):
        self.proj_state = init_lecun_normal(self.proj_state)
        nn.init.zeros_(self.proj_state.bias)
        self.proj_state2 = init_lecun_normal(self.proj_state2)
        nn.init.zeros_(self.proj_state2.bias)

    def forward(self, seq, msa, msa_full, pair, xyz_in, state, idx, use_checkpoint=False, motif_mask=None):
        """
        input:
           seq: query sequence (B, L)
           msa: seed MSA embeddings (B, N, L, d_msa)
           msa_full: extra MSA embeddings (B, N, L, d_msa_full)
           pair: initial residue pair embeddings (B, L, L, d_pair)
           xyz_in: initial BB coordinates (B, L, n_atom, 3)
           state: initial state features containing mixture of query seq, sidechain, accuracy info (B, L, d_state)
           idx: residue index
           motif_mask: bool tensor, True if motif position that is frozen, else False(L,) 
        """

        B, L = pair.shape[:2]

        if motif_mask is None:
            motif_mask = torch.zeros(L).bool()

        R_in = torch.eye(3, device=xyz_in.device).reshape(1,1,3,3).expand(B, L, -1, -1)
        T_in = xyz_in[:,:,1].clone()
        xyz_in = xyz_in - T_in.unsqueeze(-2)
        
        state = self.proj_state(state)

        R_s = list()
        T_s = list()
        alpha_s = list()
        for i_m in range(self.n_extra_block):
            R_in = R_in.detach() # detach rotation (for stability)
            T_in = T_in.detach()
            # Get current BB structure
            xyz = einsum('bnij,bnaj->bnai', R_in, xyz_in) + T_in.unsqueeze(-2)

            msa_full, pair, R_in, T_in, state, alpha = self.extra_block[i_m](msa_full, 
                                                                             pair,
                                                                             R_in, 
                                                                             T_in, 
                                                                             xyz, 
                                                                             state, 
                                                                             idx,
                                                                             motif_mask=motif_mask,
                                                                             use_checkpoint=use_checkpoint)
            R_s.append(R_in)
            T_s.append(T_in)
            alpha_s.append(alpha)

        for i_m in range(self.n_main_block):
            R_in = R_in.detach()
            T_in = T_in.detach()
            # Get current BB structure
            xyz = einsum('bnij,bnaj->bnai', R_in, xyz_in) + T_in.unsqueeze(-2)
            
            msa, pair, R_in, T_in, state, alpha = self.main_block[i_m](msa, 
                                                                       pair,
                                                                       R_in, 
                                                                       T_in, 
                                                                       xyz, 
                                                                       state, 
                                                                       idx,
                                                                       motif_mask=motif_mask,
                                                                       use_checkpoint=use_checkpoint)
            R_s.append(R_in)
            T_s.append(T_in)
            alpha_s.append(alpha)
       
        state = self.proj_state2(state)
        for i_m in range(self.n_ref_block):
            R_in = R_in.detach()
            T_in = T_in.detach()
            xyz = einsum('bnij,bnaj->bnai', R_in, xyz_in) + T_in.unsqueeze(-2)
            R_in, T_in, state, alpha = self.str_refiner(msa, 
                                                        pair, 
                                                        R_in, 
                                                        T_in, 
                                                        xyz, 
                                                        state, 
                                                        idx, 
                                                        top_k=64, 
                                                        motif_mask=motif_mask)
            R_s.append(R_in)
            T_s.append(T_in)
            alpha_s.append(alpha)

        R_s = torch.stack(R_s, dim=0)
        T_s = torch.stack(T_s, dim=0)
        alpha_s = torch.stack(alpha_s, dim=0)

        return msa, pair, R_s, T_s, alpha_s, state
