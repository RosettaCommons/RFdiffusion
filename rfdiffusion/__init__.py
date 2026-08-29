import torch

try:
    import torch_npu
    torch.npu.config.allow_internal_format = False

    # Patch torch.cdist for NPU (NPU does not support cdist natively)
    _orig_cdist = torch.cdist
    def _npu_cdist(x1, x2, p=2.0, compute_mode='use_mm_for_euclid_dist_if_necessary', **kwargs):
        if x1.device.type == 'npu' or (x2 is not None and hasattr(x2, 'device') and x2.device.type == 'npu'):
            if x1.dim() == 2:
                x1 = x1.unsqueeze(0)
                x2 = x2.unsqueeze(0)
                squeeze = True
            else:
                squeeze = False
            x1_sq = (x1 * x1).sum(dim=-1, keepdim=True)
            x2_sq = (x2 * x2).sum(dim=-1, keepdim=True)
            cross = torch.bmm(x1, x2.transpose(-1, -2))
            dist_sq = x1_sq + x2_sq.transpose(-1, -2) - 2 * cross
            dist_sq = dist_sq.clamp(min=0)
            result = torch.sqrt(dist_sq + 1e-12)
            if squeeze:
                result = result.squeeze(0)
            return result
        return _orig_cdist(x1, x2, p=p, compute_mode=compute_mode, **kwargs)
    torch.cdist = _npu_cdist
except ImportError:
    pass
