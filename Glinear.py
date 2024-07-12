import numpy as np
import scipy.sparse as sp

def sys_matrix(nx, ny, nb, na, mask):
#function G = Glinear(nx, ny, nb, na, mask)
#
#	Generate a geometric system matrix for tomographic projection
#	based on simple linear interpolation.
#	Contains exactly two g_{ij}'s per pixel per projection angle.
#	This system model is pretty inadequate for reconstructing
#	real tomographic data, but is useful for simple simulations.
#
#	The image size is nx * ny.
#	The sinogram size is nb * na (n_radial_bins X n_angles).
#	Returned G is [nb * na, nx * ny] - this size is need for wtfmex saving.
#
#	Copyright Apr 2000, Jeff Fessler, University of Michigan

#if nargin < 1, nx = 32; end
#if nargin < 2, ny = nx; end
#if nargin < 3, nb = nx; end
#if nargin < 4, na = floor(nb * pi/2); end
#if nargin < 5, mask = logical(ones(nx,ny)); end

#
#	default: run a test demo
#
#if nargin < 1
	#x = shepplogan(nx, ny, 1);
	#ix = [-(nx-1)/2:(nx-1)/2]';
	#iy = [-(ny-1)/2:(ny-1)/2];
	#rr = sqrt(outer_sum(ix.^2, iy.^2));
	#mask = rr < nx/2-1;
	#clf
	#G = Glinear(nx, ny, nb, na, mask);
	#y = G * x(:);		% forward projection
	#y = reshape(y, nb, na);	% reshape into sinogram array
	#sino = zeros(nb, na);	sino(nb/2, 10) = 1;
	#b = reshape(G' * sino(:), nx, ny);
	#im(221, mask, 'support mask')
	#im(222, x, 'test image')
	#im(223, y, 'sinogram'), xlabel ib, ylabel ia
	#im(224, b, 'backproject 1 ray')
#return
#end

#
#	pixel centers
#
    x,y = np.meshgrid(np.arange(0,nx)-(nx-1)/2., np.arange(0,ny)-(ny-1)/2.)
    x = x[mask]
    y = y[mask]
    npts = len(x)		# sum(mask(:)) - total # of support pixels

    angle = np.arange(0,na) * (np.pi/na)
# [na,np] projected pixel center
    tau = np.outer(np.cos(angle),x) + np.outer(np.sin(angle), y)
    tau += nb/2.		# counting from 0 (python)
    ibl = np.floor(tau).astype(int)		# left bin
    val = 1 - (tau-ibl)		# weight value for left bin
    #print(val)
#ii is the absolute sinogram index (from 0 to na*np) in image space
    ii = ibl + nb*np.outer(np.arange(0,na),np.ones(npts)).astype(int)

    #good = ones ([na,npts])
    good = np.logical_and(ibl > -1,ibl < nb)	# within FOV cases
    if not(good.any()):
        print ('FOV too small')
        return

    #print(val[good])
#np = sum(mask(:));
#nc = np;	jj = 1:np;
# compact G
    nc = nx * ny
    jj = mask.ravel().nonzero()
# all-column G
    jj=np.outer(np.ones(na),jj).astype(int)

#left bin
    G1 = sp.csc_matrix((val[good],(ii[good],jj[good])),shape=[nb*na,nc])
#right bin

#redefine valid indices
    good1 = (ibl > -2) & (ibl < nb-1)	# within FOV cases, but shifted by one
    G2 = sp.csc_matrix((1-val[good1],(ii[good1]+1,jj[good1])),shape=[nb*na,nc])

    return G1 + G2