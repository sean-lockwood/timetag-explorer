"""STIS FUV geometric distortion correction: a fourth-order polynomial in the detector position.

THE COEFFICIENTS BELOW ARE PLACEHOLDERS (all zero), so `gdcorr` returns its input unchanged. While they are zero
`PLACEHOLDER` is True, TIME-TAG Explorer leaves the geometric-distortion switch off and disabled, and nothing is
corrected. Real coefficients, taken from a public reference file, will replace them in a later release.
"""
a = [0.0] * 14
b = [0.0] * 14

# True while the coefficients are placeholders; becomes False by itself once any of them is non-zero.
PLACEHOLDER = not (any(a) or any(b))


def polyd4(x, y, c):
   ff = c[0]*y       + c[1]*x        + \
        c[2]*y*y     + c[3]*x*y      + c[4]*x*x       + \
        c[5]*y*y*y   + c[6]*y*y*x    + c[7]*y*x*x     + c[8]*x*x*x    + \
        c[9]*y*y*y*y + c[10]*y*y*y*x + c[11]*y*y*x*x  + c[12]*y*x*x*x + c[13]*x*x*x*x
   return ff


def gdcorr(xr, yr):
    '''
    Apply the STIS FUV geometric distortion correction.

    INPUTS:
        xr, yr   - Input raw (x, y) positions

    OUTPUTS:
        xc, yc   - Output corrected (x, y) positions
    '''
    ndx, ndy = 512., 512.

    xrn = xr/2. - ndx
    yrn = yr/2. - ndy
    ndx = polyd4(xrn, yrn, a)
    ndy = polyd4(xrn, yrn, b)
    xc  = xr + ndx*2.
    yc  = yr + ndy*2.

    return xc, yc
