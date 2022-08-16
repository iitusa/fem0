package math0;

import org.apache.commons.math.linear.*;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

@RunWith(Parameterized.class)
public class TestFEF0 {

    static final private int nnn = 720;
    static final double __xu = 7.0;
    static final double __xl = - __xu;
    static final double theEPS = 0.01;

    @Parameterized.Parameters
    public static List<Object[]> parameters() {
        final List<Object[]> params = new ArrayList<>();
        for (int kk=0; kk<=10; ++kk) {
            final double ccc = -2.0 + 0.4 * kk;
            for (int ii=0; ii<10; ++ii) {
                final double aaa = -1.8 + 0.4 * ii;
                if (ccc>0) {
                    final Object[] sub = new Object[3];
                    sub[0] = aaa;
                    sub[1] = Double.NaN;
                    sub[2] = ccc;
                    params.add(sub);
                } else {
                    for (int jj=0; jj<10; ++jj) {
                        final double bbb = -1.8 + 0.4 * jj;
                        final Object[] sub = new Object[3];
                        sub[0] = aaa;
                        sub[1] = bbb;
                        sub[2] = ccc;
                        params.add(sub);
                    }
                }
            }
        }
        return params;
    }

    final private Fem fem = new Fem(nnn);
//    private RealMatrix qqq;
    final double aaa;
    final double bbb;
    final double ccc;
//    private RealVector rv;
//    final private double[] expected = new double[nnn+1];


    public TestFEF0(double aaa, double bbb, double ccc) {
        this.aaa = aaa;
        this.bbb = bbb;
        this.ccc = ccc;
    }

//    @Before
//    public void setUp() {
//        qqq = fem.makeTri(nnn-1, __xl, __xu);
//    }

    private RealVector setupExpected(int nnn, double xl, double xu, double aaa, double bbb, double ccc) {
        final double[] expected = new double[nnn +1];
        for (int ii=0; ii< nnn +1; ++ii) {
            final double xxx = ((nnn - ii) * xl + ii * xu) / nnn;
            if (ccc>0) {
                expected[ii] = aaa * Math.exp(Math.sqrt(ccc) * xxx);
            } else if (ccc<0) {
                expected[ii] = aaa * Math.cos(Math.sqrt(-ccc) * xxx) + bbb * Math.sin(Math.sqrt(-ccc) * xxx);
            } else {
                expected[ii] = aaa * xxx + bbb;
            }
        }
        return MatrixUtils.createRealVector(expected);
    }

//    @Test
//    @Ignore
//    public void test9() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
//        final EigenDecomposition eigenDecomposition = new EigenDecompositionImpl(qqq.getSubMatrix(1,nnn-2,1,nnn-2), Double.MIN_NORMAL);
//        Logger.getAnonymousLogger().info(Arrays.toString(eigenDecomposition.getRealEigenvalues()));
//    }

//    @Test
//    @Ignore
//    public void printQQQ1() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
//        final RealMatrix qqqqq = MatrixUtils.createRealIdentityMatrix(nnn-1).scalarMultiply(ccc).add(qqq);
//        Logger.getAnonymousLogger().info("" + qqqqq);
//    }

    @Test
    public void test1b() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
        final RealVector expected = setupExpected(nnn, __xl, __xu, aaa, bbb, ccc);
        final RealVector fef = fem.solve(__xl, __xu, expected.getEntry(0), expected.getEntry(nnn), ccc);

        Logger.getAnonymousLogger().fine("EXPECT" + expected);
        Logger.getAnonymousLogger().fine("ACTUAL" + fef);
        for (int ii=0; ii<fef.getDimension(); ++ii) {
            final double actual = fef.getEntry(ii);
            final double expect = expected.getEntry(ii+1);
            final double eps = Math.max(theEPS * Math.sqrt(actual*actual + expect*expect), theEPS);
            Assert.assertEquals("A="+ aaa + "\tB="+ bbb + "\tC="+ ccc + "\tii=" + ii, expect, actual, eps);
        }
    }

    @Test
    public void test1a() throws NotSymmetricMatrixException, NotPositiveDefiniteMatrixException {
        final double delta = (__xu-__xl) / nnn;
        final double delta2 = delta*delta;

        final RealVector expected = setupExpected(nnn, __xl, __xu, aaa, bbb, ccc);

        final double[] actuals = new double[nnn-1];
        for (int ii=0; ii<nnn-1; ++ii) {
            actuals[ii] = (expected.getEntry(ii) + expected.getEntry(ii+2) - 2.0 * expected.getEntry(ii+1)) / delta2;
        }
        final RealVector realVector2 = MatrixUtils.createRealVector(actuals);
        Logger.getAnonymousLogger().fine("" + realVector2);
        for (int ii=0; ii<actuals.length; ++ii) {
            final double expect = ccc * expected.getEntry(ii+1);
            final double actual = actuals[ii];
            final double eps = Math.max(theEPS * Math.sqrt(actual*actual + expect*expect), theEPS);
            Assert.assertEquals(""+ ii, expect, actual, eps);
        }
    }

//    @Test
//    @Ignore
//    public void test2() {
//        final double x0 = Math.PI * 0.75;
//        final double eps=1e-6;
//        final double y0 = Math.cos(x0);
//        final double y1 = Math.cos(x0-eps);
//        final double y2 = Math.cos(x0 + eps);
//        Logger.getAnonymousLogger().info("" + (y1+y2-2.0*y0)/eps/eps);
//    }
//
//    @Test
//    @Ignore
//    public void test2a() {
//        final double x0 = Math.PI * 0.75;
//        final double eps=1e-6;
//        final double y0 = Math.cos(x0);
//        final double y1 = Math.cos(x0-eps);
//        final double y2 = Math.cos(x0+eps);
//        Logger.getAnonymousLogger().info("" + (2.0 * y2-2.0*y0)/eps/eps);
//    }
}
