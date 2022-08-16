package math0;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.random.GaussianRandomGenerator;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.logging.Logger;

public class TestQR0 {

    static private final int nnn = 199;

    static private JDKRandomGenerator jrg;
    static private GaussianRandomGenerator grg;
    static private RealMatrix mmm;
    static private RealMatrix pdm;

    @BeforeClass
    static public void setUp() {

        jrg = new JDKRandomGenerator();
        jrg.setSeed(System.currentTimeMillis());
        grg = new GaussianRandomGenerator(jrg);

        final double[][] vvv = new double[nnn][nnn];
        for (int ii=0; ii<nnn; ++ii) for (int jj=0; jj<nnn; ++jj) vvv[ii][jj] = grg.nextNormalizedDouble() / Math.sqrt( nnn / Math.PI);
        mmm = MatrixUtils.createRealMatrix(vvv);
        Logger.getAnonymousLogger().fine("" + mmm);
        pdm = mmm.transpose().multiply(mmm);
    }

    @Test
    public void testCholesky() {
        final CholeskyDecomposition cholesky = new CholeskyDecomposition(pdm);
        Logger.getAnonymousLogger().info("" + cholesky.getDeterminant());
    }

    @Test
    public void testQR() {

        final EigenDecomposition eigenDecomposition = new EigenDecomposition(pdm, Double.MIN_NORMAL);
        Logger.getAnonymousLogger().info(Arrays.toString(eigenDecomposition.getRealEigenvalues()));
        Logger.getAnonymousLogger().info("" + eigenDecomposition.getDeterminant());

        final RealMatrix unitary = eigenDecomposition.getV().multiply(eigenDecomposition.getVT());
        for (int ii=0; ii<nnn; ++ii) for (int jj=0; jj<nnn; ++jj) Assert.assertEquals("" + ii + "-" + jj, ii==jj ? 1.0 : 0.0, unitary.getEntry(ii, jj), 1e-6);
    }

    @Test
    public void testQR0() {
        final double[][] doubles = {{1.0,4.0},{4.0,16.0}};
        final RealMatrix nnn = MatrixUtils.createRealMatrix(doubles);
        final QRDecomposition qr = new QRDecomposition(nnn);
        Logger.getAnonymousLogger().info("H=" + qr.getH());
        Logger.getAnonymousLogger().info("Q=" + qr.getQ());
        Logger.getAnonymousLogger().info("R=" + qr.getR());
    }
}
