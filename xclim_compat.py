"""
xclim_compat.py  -  restore the one xarray internal that xclim 0.54.0 uses and
xarray 2026.7.0 no longer provides.

The failure. xclim/indices/generic.py:get_op ends with

    return xr.core.ops.get_op(binary_op)

and xarray 2026.7.0 has no xarray.core.ops module at all; importing it raises
ModuleNotFoundError, and the attribute access raises AttributeError. Anything
that compares an array to a threshold goes through that call, so percentile_doy
fails immediately and every percentile index with it.

xclim maps symbols to names first, through
xclim.indices.generic.binary_ops:

    {">": "gt", "<": "lt", ">=": "ge", "<=": "le", "==": "eq", "!=": "ne"}

so the string reaching xarray is always one of gt, lt, ge, le, eq, ne. Those are
exactly the names in the operator module, and xarray's old get_op returned the
operator function for the name. The shim below does the same thing for that
fixed set and raises on anything else, rather than forwarding arbitrary
attribute lookups.

The point of pinning xclim to 0.54.0 is
that the hub runs 0.54.0, so the Derecho numbers can be compared with the bucket
numbers without a library difference in the way. Upgrading xclim here would
remove that guarantee for the sake of a two-line incompatibility. The cleaner
fix is to pin xarray to whatever the hub runs, and that is the fix I would
prefer: check xr.__version__ on the hub and pin the Casper environment to match.
Until that is known, this keeps the two xclim versions identical and changes
nothing about what xclim computes.

The shim prints when it patches, and it verifies the
six operators against the operator module before declaring success.

"""

import operator
import sys
import types

_NAMES = ('gt', 'lt', 'ge', 'le', 'eq', 'ne')


def _get_op(name):
    if name not in _NAMES:
        raise ValueError(
            f'xclim_compat shim: unexpected operator {name!r}. The shim covers '
            f'{_NAMES}, which is what xclim.indices.generic.binary_ops maps to. '
            f'A name outside that set means xclim is using xarray.core.ops for '
            f'something this shim was not written for; pin xarray to a version '
            f'that still provides it instead of extending this.')
    return getattr(operator, name)


def patch(verbose=True):
    """
    Install xarray.core.ops.get_op if xarray no longer provides it.

    Returns True if a patch was applied, False if xarray already has it.
    
    """
    import xarray as xr
    import xarray.core as xrcore

    if hasattr(xrcore, 'ops') and hasattr(xrcore.ops, 'get_op'):
        return False

    mod = types.ModuleType('xarray.core.ops')
    mod.get_op = _get_op
    sys.modules['xarray.core.ops'] = mod
    xrcore.ops = mod

    # Verify before trusting it: these are the comparisons every percentile
    # index depends on, and a wrong one would produce plausible, wrong counts.
    checks = [('gt', 2, 1, True), ('lt', 2, 1, False), ('ge', 1, 1, True),
              ('le', 1, 2, True), ('eq', 1, 1, True), ('ne', 1, 1, False)]
    for name, a, b, want in checks:
        got = xr.core.ops.get_op(name)(a, b)
        if got is not want:
            raise RuntimeError(f'xclim_compat shim failed self-test on {name}: '
                               f'{a} {name} {b} gave {got}, expected {want}')

    if verbose:
        print(f'  xclim_compat: patched xarray.core.ops.get_op for '
              f'xarray {xr.__version__} (xclim 0.54.0 expects the removed '
              f'module). Self-test passed on {", ".join(_NAMES)}.')
    return True
