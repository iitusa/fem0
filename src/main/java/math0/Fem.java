package math0;

import org.apache.commons.math3.linear.*;

/**
 * finite element method
 */
public class Fem {

    final private int nnn;

    public Fem(int nnn) {
        this.nnn = nnn;
    }

    /**
     * create a band diagonal matrix for d^2u=u differential equation
     * @param nnn size
     * @param xl lower bound for calculating delta
     * @param xu upper bound for calculating delta
     * @return RealMatrix with a band diagonal
     */
    public RealMatrix makeTri(int nnn, double xl, double xu) {
        final double delta = (xu-xl) / (nnn+1);
        final double delta2 = 1.0/delta/delta;
        final double[][] mmm = new double[nnn][nnn];
        for (int ii=0; ii<mmm.length; ++ii) for (int jj=0; jj<mmm.length; ++jj) {
            if (ii==jj) {
                mmm[ii][jj] = - 2.0 * delta2;
            } else if (Math.abs(ii-jj) == 1) {
                mmm[ii][jj] = delta2;
            } else {
                mmm[ii][jj] = 0.0;
            }
        }
        return MatrixUtils.createRealMatrix(mmm);
    }

    /**
     * solve the differential equation
     * @param xl
     * @param xu
     * @param yl
     * @param yu
     * @param ccc
     * @return
     * @throws NotSymmetricMatrixException
     * @throws NotPositiveDefiniteMatrixException
     */
    public RealVector solve(double xl, double xu, double yl, double yu, double ccc) {

        final RealMatrix qqq = this.makeTri(nnn - 1, xl, xu);

        final double delta = (xu-xl) / nnn;
        final double delta2 = 1.0/delta/delta;

        final double[] vvv = new double[nnn-1];
        for (int kk=0; kk<vvv.length; ++kk) {
            if (kk==0) {
                vvv[kk] = yl * delta2;
            } else if (kk==nnn-2) {
                vvv[kk] = yu * delta2;
            } else {
                vvv[kk] = 0.0;
            }
        }
        final RealVector rv = MatrixUtils.createRealVector(vvv);

        final RealMatrix theQQQ = MatrixUtils.createRealIdentityMatrix(nnn-1).scalarMultiply(ccc).subtract(qqq);
        final RealVector fef;
        if (ccc>=0) {
            final CholeskyDecomposition lu = new CholeskyDecomposition(theQQQ);
            fef = lu.getSolver().solve(rv);
        } else {
            final LUDecomposition lu = new LUDecomposition(theQQQ);
            fef = lu.getSolver().solve(rv);
        }
        return fef;
    }
}
