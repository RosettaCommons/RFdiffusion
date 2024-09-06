import warnings

try:
    from torch._C import _nvtx
    from torch.cuda.nvtx import range as nvtx_range
except ImportError as e:
    warnings.warn(f'NVTX is not available: {e}')

    class MockNvtx:
        def __enter__(self):
            pass
        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    nvtx_range = lambda t: MockNvtx()
