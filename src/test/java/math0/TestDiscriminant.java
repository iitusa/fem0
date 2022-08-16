package math0;

import org.apache.commons.math.linear.*;
import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math.stat.descriptive.moment.FirstMoment;
import org.apache.commons.math.stat.regression.OLSMultipleLinearRegression;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.logging.Logger;

public class TestDiscriminant {

    static private final double[][] zz0 = {
            {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
            {178, 169.8, 165.1, 159.4, 163,   160.1, 178.4, 177.5, 172.4, 172.1,
             180, 171.2, 173.0, 179.0, 170.3, 174,   174.3, 168.7, 171.2, 170.6},
            {56, 52, 56, 44, 54, 51.5, 60.5, 62, 64.5, 65, 78, 58, 67, 78, 55, 65, 71, 59, 72.5, 56},
            {80, 87.5, 82.5, 73, 83, 80, 83, 85, 87, 85, 97, 84.5, 92.8, 98, 85.6, 89, 88.5, 86.6, 97, 85},
            {93.3, 88,   87.5, 85.2, 88.5, 86.6, 91.7, 92.4, 92.4, 91.7,
             98.5, 91.5, 91.9, 94.8, 88.5, 93.5, 91.5, 93,   92.2, 91}
    };

    static private RealMatrix xx0;
    static private RealMatrix xx0a;
    static private RealMatrix xx1;
    static private RealMatrix xx2;

    @Before
    public void setup() {
        xx0a = MatrixUtils.createRealMatrix(zz0).transpose();
        xx0 = xx0a.getSubMatrix(0,19,1,4);
        xx1 = xx0.getSubMatrix(0, 9, 0, 3);
        xx2 = xx0.getSubMatrix(10, 19, 0, 3);
    }

    static private volatile double _vvv;

    @Test
    @Ignore
    public void cov1() {
        final Covariance cov1 = new Covariance(xx1);
        Logger.getAnonymousLogger().info("" + cov1.getCovarianceMatrix());
    }

    @Test
    @Ignore
    public void cov2() {
        final Covariance cov2 = new Covariance(xx2);
        Logger.getAnonymousLogger().info("" + cov2.getCovarianceMatrix());
    }

    @Test
    @Ignore
    public void cov0() {
        Logger.getAnonymousLogger().fine("" + xx0);
        final Covariance cov0 = new Covariance(xx0);
        Logger.getAnonymousLogger().info("" + cov0.getCovarianceMatrix());
    }

    @Test
    public void mahalanobis() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {

        final FirstMoment fm = new FirstMoment();
        final double[] diff = new double[xx0.getColumnDimension()];
        for (int ii=0; ii<diff.length; ++ii) {
            fm.setData(xx1.getColumn(ii));
            fm.evaluate();
            final double avg1 = fm.getResult();
            fm.setData(xx2.getColumn(ii));
            fm.evaluate();
            final double avg2 = fm.getResult();
            diff[ii] = avg2 - avg1;
        }
        Logger.getAnonymousLogger().info(Arrays.toString(diff));

        RealMatrix cov0 = new Covariance(xx0).getCovarianceMatrix();
        {
        }
        final CholeskyDecomposition cholesky0 = new CholeskyDecompositionImpl(cov0);
        final double[] solve = cholesky0.getSolver().solve(diff);
        final double vvv = new ArrayRealVector(diff).dotProduct(solve);
        Logger.getAnonymousLogger().info("" + vvv + "\t" + vvv*vvv + "\t" + Math.sqrt(vvv));
        _vvv = vvv*vvv;
    }

    @Test
    public void fit0() {
        final double sqrt = 3.144285975645797; // Math.sqrt(11.5);
        final double[] yy0 = new double[20];
        for (int ii=0; ii<yy0.length; ++ii) {
            yy0[ii] = 2*ii<yy0.length ? sqrt : -sqrt;
        }

        final OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
        ols.setNoIntercept(false);
        ols.newSampleData(yy0, xx0.getData());
        final double[] doubles = ols.estimateRegressionParameters();
        Logger.getAnonymousLogger().info(Arrays.toString(doubles));
        final double[] doubles2 = ols.estimateResiduals();
        for (int ii=0; ii<doubles2.length; ++ii) {
            doubles2[ii] = yy0[ii] - doubles2[ii];
            doubles2[ii] *= 1000;
            doubles2[ii] = Math.floor(doubles2[ii]);
            doubles2[ii] /= 1000;
        }
        Logger.getAnonymousLogger().info(Arrays.toString(doubles2));
    }

    @Test
    public void fit1() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
        final double sqrt = 3.144285975645797; // Math.sqrt(11.5);
        final double[] yy0 = new double[20];
        for (int ii=0; ii<yy0.length; ++ii) {
            yy0[ii] = 2*ii<yy0.length ? sqrt : -sqrt;
        }

        final RealMatrix cov0 = xx0a.transpose().multiply(xx0a);
        final CholeskyDecomposition cholesky0 = new CholeskyDecompositionImpl(cov0);
        final double[][] xxx0 = new double[1][];
        xxx0[0] = yy0;
        final RealMatrix solve = cholesky0.getSolver().solve(MatrixUtils.createRealMatrix(xxx0).multiply(xx0a).transpose());
        Logger.getAnonymousLogger().info("" + solve.transpose());
    }

    @Test
    public void fit2() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
        final double sqrt = 0.19;
        final double[] yy0 = new double[20];
        for (int ii=0; ii<yy0.length; ++ii) {
            yy0[ii] = 2*ii<yy0.length ? sqrt : -sqrt;
//            yy0[ii] += 47.349;
//            yy0[ii] += -10;
        }

        final RealMatrix xx1a = xx0a.getSubMatrix(0,9,0,4);
        final RealMatrix xx2a = xx0a.getSubMatrix(10, 19, 0, 4);
        final RealMatrix cov1 = xx1a.transpose().multiply(xx1a);
        final RealMatrix cov2 = xx2a.transpose().multiply(xx2a);
        final RealMatrix cov0 = cov1.add(cov2).scalarMultiply(1.0/18.0);
        final CholeskyDecomposition cholesky0 = new CholeskyDecompositionImpl(cov0);
        final double[][] xxx0 = new double[1][];
        xxx0[0] = yy0;
        final RealMatrix solve = cholesky0.getSolver().solve(MatrixUtils.createRealMatrix(xxx0).multiply(xx0a).transpose());
        Logger.getAnonymousLogger().info("" + solve.transpose());

//        final RealVector realVector = MatrixUtils.createRealVector(yy0);
//        final RealMatrix solve2 = cholesky0.getSolver().solve(realVector.createRealMatrix(xxx0).multiply(xx0a).transpose());
    }
}
