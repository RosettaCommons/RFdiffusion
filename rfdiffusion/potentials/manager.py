import torch
from rfdiffusion.potentials import potentials as potentials
import numpy as np 


def make_contact_matrix(nchain, intra_all=False, inter_all=False, contact_string=None):
    """
    Calculate a matrix of inter/intra chain contact indicators
    
    Parameters:
        nchain (int, required): How many chains are in this design 
        
        contact_str (str, optional): String denoting how to define contacts, comma delimited between pairs of chains
            '!' denotes repulsive, '&' denotes attractive
    """
    alphabet   = [a for a in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
    letter2num = {a:i for i,a in enumerate(alphabet)}
    
    contacts   = np.zeros((nchain,nchain))
    written    = np.zeros((nchain,nchain))
    
    
    # intra_all - everything on the diagonal has contact potential
    if intra_all:
        contacts[np.arange(nchain),np.arange(nchain)] = 1
    
    # inter all - everything off the diagonal has contact potential
    if inter_all:
        mask2d = np.full_like(contacts,False)
        for i in range(len(contacts)):
            for j in range(len(contacts)):
                if i!=j:
                    mask2d[i,j] = True
        
        contacts[mask2d.astype(bool)] = 1


    # custom contacts/repulsions from user 
    if contact_string != None:
        contact_list = contact_string.split(',') 
        for c in contact_list:
            assert len(c) == 3
            i,j = letter2num[c[0]],letter2num[c[2]]

            symbol = c[1]

            assert symbol in ['!','&']
            if symbol == '!':
                contacts[i,j] = -1
                contacts[j,i] = -1
            else:
                contacts[i,j] = 1
                contacts[j,i] = 1
            
    return contacts 


def calc_nchains(symbol, components=1):
    """
    Calculates number of chains for given symmetry 
    """
    S = symbol.lower()

    if S.startswith('c'):
        return int(S[1:])*components 
    elif S.startswith('d'):
        return 2*int(S[1:])*components 
    elif S.startswith('o'):
        raise NotImplementedError()
    elif S.startswith('t'):
        return 12*components
    else:
        raise RuntimeError('Unknown symmetry symbol ',S)


class PotentialManager:
    '''
        Class to define a set of potentials from the given config object and to apply all of the specified potentials
        during each cycle of the inference loop.

        Author: NRB
    '''

    def __init__(self, 
                 potentials_config, 
                 ppi_config, 
                 diffuser_config, 
                 inference_config,
                 hotspot_0idx,
                 binderlen,
                 ):

        self.potentials_config = potentials_config
        self.ppi_config        = ppi_config
        self.inference_config  = inference_config

        self.guide_scale = potentials_config.guide_scale
        self.guide_decay = potentials_config.guide_decay
    
        if potentials_config.guiding_potentials is None: 
            setting_list = []
        else: 
            setting_list = [self.parse_potential_string(potstr) for potstr in potentials_config.guiding_potentials]


        # PPI potentials require knowledge about the binderlen which may be detected at runtime
        # This is a mechanism to still allow this info to be used in potentials - NRB 
        if binderlen > 0:
            binderlen_update   = { 'binderlen': binderlen }
            hotspot_res_update = { 'hotspot_res': hotspot_0idx }

            for setting in setting_list:
                if setting['type'] in potentials.require_binderlen:
                    setting.update(binderlen_update)

        self.potentials_to_apply = self.initialize_all_potentials(setting_list)
        self.T = diffuser_config.T
        
    def is_empty(self):
        '''
            Check whether this instance of PotentialManager actually contains any potentials
        '''

        return len(self.potentials_to_apply) == 0

    def parse_potential_string(self, potstr):
        '''
            Parse a single entry in the list of potentials to be run to a dictionary of settings for that potential.

            An example of how this parsing is done:
            'setting1:val1,setting2:val2,setting3:val3' -> {setting1:val1,setting2:val2,setting3:val3}
        '''

        setting_dict = {entry.split(':')[0]:entry.split(':')[1] for entry in potstr.split(',')}

        for key in setting_dict:
            if not key == 'type': setting_dict[key] = float(setting_dict[key])

        return setting_dict

    def initialize_all_potentials(self, setting_list):
        '''
            Given a list of potential dictionaries where each dictionary defines the configurations for a single potential,
            initialize all potentials and add to the list of potentials to be applies
        '''

        to_apply = []

        for potential_dict in setting_list:
            assert(potential_dict['type'] in potentials.implemented_potentials), f'potential with name: {potential_dict["type"]} is not one of the implemented potentials: {potentials.implemented_potentials.keys()}'

            kwargs = {k: potential_dict[k] for k in potential_dict.keys() - {'type'}}

            # symmetric oligomer contact potential args
            if self.inference_config.symmetry:

                num_chains = calc_nchains(symbol=self.inference_config.symmetry, components=1) # hard code 1 for now 
                contact_kwargs={'nchain':num_chains,
                                'intra_all':self.potentials_config.olig_intra_all,
                                'inter_all':self.potentials_config.olig_inter_all,
                                'contact_string':self.potentials_config.olig_custom_contact }
                contact_matrix = make_contact_matrix(**contact_kwargs)
                kwargs.update({'contact_matrix':contact_matrix})


            to_apply.append(potentials.implemented_potentials[potential_dict['type']](**kwargs))

        return to_apply

    def compute_all_potentials(self, xyz):
        '''
            This is the money call. Take the current sequence and structure information and get the sum of all of the potentials that are being used
        '''

        potential_list = [potential.compute(xyz) for potential in self.potentials_to_apply]
        potential_stack = torch.stack(potential_list, dim=0)

        return torch.sum(potential_stack, dim=0)

    def get_guide_scale(self, t):
        '''
        Given a timestep and a decay type, get the appropriate scale factor to use for applying guiding potentials
        
        Inputs:
        
            t (int, required):          The current timestep
        
        Output:
        
            scale (int):                The scale factor to use for applying guiding potentials
        
        '''
        
        implemented_decay_types = {
                'constant': lambda t: self.guide_scale,
                # Linear interpolation with y2: 0, y1: guide_scale, x2: 0, x1: T, x: t
                'linear'  : lambda t: t/self.T * self.guide_scale,
                'quadratic' : lambda t: t**2/self.T**2 * self.guide_scale,
                'cubic' : lambda t: t**3/self.T**3 * self.guide_scale
        }
        
        if self.guide_decay not in implemented_decay_types:
            sys.exit(f'decay_type must be one of {implemented_decay_types.keys()}. Received decay_type={self.guide_decay}. Exiting.')
        
        return implemented_decay_types[self.guide_decay](t)


        
