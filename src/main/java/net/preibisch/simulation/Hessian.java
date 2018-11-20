package net.preibisch.simulation;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import net.imglib2.Cursor;
import net.imglib2.Positionable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Sampler;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

public class Hessian
{
	/**
	 * @param randomAccessible
	 * @return - the Math.abs( largest ) eigenvalue as n-dim image (A), and the corresponding eigenvector as (n+1)-dim image (B)
	 */
	public static Pair< Img< FloatType >, Img< FloatType > > largestEigenVector( final RandomAccessibleInterval< FloatType > randomAccessible )
	{
		final int n = randomAccessible.numDimensions();

		// blur a bit
		Gauss3.gauss( 1.0, Views.extendMirrorSingle( randomAccessible ), randomAccessible );

		final long[] dimEigenValue = new long[ n ];
		randomAccessible.dimensions( dimEigenValue );

		final long[] dimEigenVector = new long[ n + 1 ];

		for ( int d = 0; d < n; ++d )
			dimEigenVector[ d ] = dimEigenValue[ d ];

		dimEigenVector[ n ] = 3; // vector

		final Img< FloatType > eigenVec = new ArrayImgFactory<>( new FloatType() ).create( dimEigenVector );
		final Img< FloatType > eigenVal = new ArrayImgFactory<>( new FloatType() ).create( dimEigenValue );

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		final Cursor< FloatType > cursor = Views.iterable( randomAccessible ).localizingCursor();
		final RandomAccess< FloatType > raEVal = eigenVal.randomAccess();
		final RandomAccess< FloatType > raEVec = eigenVec.randomAccess();
		final long[] posN = new long[ eigenVec.numDimensions() ];

		final RandomAccess< FloatType > ra = Views.extendMirrorSingle( randomAccessible ).randomAccess();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( posN );
			posN[ n ] = 0; // n+1 dim

			ra.setPosition( cursor );

			computeHessianMatrix3D( ra, matrix );

			final double eigenValue = computeLargestEigenVectorAndValue3d( matrix, eigenVector );

			raEVec.setPosition( posN );
			raEVec.get().set( (float)eigenVector[ 0 ] );

			raEVec.fwd( n );
			raEVec.get().set( (float)eigenVector[ 1 ] );

			raEVec.fwd( n );
			raEVec.get().set( (float)eigenVector[ 2 ] );

			raEVal.setPosition( cursor );
			raEVal.get().set( (float)eigenValue );
		}

		return new ValuePair< Img<FloatType>, Img<FloatType> >( eigenVal, eigenVec );
	}

	/**
	 * 
	 * @param hessianMatrix
	 * @param eigenVector - fills up the eigenvector
	 * @return eigenvalue
	 */
	public static double computeLargestEigenVectorAndValue3d( final double[][] hessianMatrix, final double[] eigenVector )
	{
		final Matrix m = new Matrix( hessianMatrix );
		final EigenvalueDecomposition evDec = new EigenvalueDecomposition( m );

		final double[] result = evDec.getImagEigenvalues();

		for ( int i = 0; i < result.length; ++i )
			if ( result[i] > 0 )
			{
				// found imaginary parts
				eigenVector[ 0 ] = 1;
				eigenVector[ 1 ] = eigenVector[ 2 ] = 0;
				return 0; 
			}

		final double[] ev = evDec.getRealEigenvalues();

		int biggestEigenValueIndex = 0;
		double biggestEigenValue = ev[ 0 ];

		for ( int i = 1; i < ev.length; ++i )
		{
			if ( Math.abs( ev[ i ] ) > Math.abs( biggestEigenValue ) )
			{
				biggestEigenValue = ev[ i ];
				biggestEigenValueIndex = i;
			}
		}

		final Matrix evec = evDec.getV();

		eigenVector[ 0 ] = evec.get( 0, biggestEigenValueIndex );
		eigenVector[ 1 ] = evec.get( 1, biggestEigenValueIndex );
		eigenVector[ 2 ] = evec.get( 2, biggestEigenValueIndex );

		return biggestEigenValue;
	}

	/**
	 * computes the Hessian Matrix in local 3x3x3 environment. Works with RandomAccess and RealRandomAccess
	 * 
	 * @param ra
	 * @param hessianMatrix
	 */
	public static < A extends Positionable & Sampler< FloatType > > void computeHessianMatrix3D( final A ra, final double[][] hessianMatrix )
	{
		final double temp = 2 * ra.get().get();

		// xx
		//hessianMatrix[0][0] = img.get(x + 1, y, z) - temp + img.get(x - 1, y, z);
		ra.fwd( 0 );
		hessianMatrix[0][0] = ra.get().get();
		hessianMatrix[0][0] -= temp;
		ra.bck( 0 );
		ra.bck( 0 );
		hessianMatrix[0][0] += ra.get().get();
		ra.fwd( 0 );

		// yy
		//hessianMatrix[1][1] = img.get(x, y + 1, z) - temp + img.get(x, y - 1, z);
		ra.fwd( 1 );
		hessianMatrix[1][1] = ra.get().get();
		hessianMatrix[1][1] -= temp;
		ra.bck( 1 );
		ra.bck( 1 );
		hessianMatrix[1][1] += ra.get().get();
		ra.fwd( 1 );

		// zz
		//hessianMatrix[2][2] = img.get(x, y, z + 1) - temp + img.get(x, y, z - 1);
		ra.fwd( 2 );
		hessianMatrix[2][2] = ra.get().get();
		hessianMatrix[2][2] -= temp;
		ra.bck( 2 );
		ra.bck( 2 );
		hessianMatrix[2][2] += ra.get().get();
		ra.fwd( 2 );

		// xy
		ra.fwd( 0 );
		ra.fwd( 1 );
		double a = ra.get().get();
		ra.bck( 0 );
		ra.bck( 0 );
		double b = ra.get().get();
		ra.fwd( 0 );
		ra.fwd( 0 );
		ra.bck( 1 );
		ra.bck( 1 );
		double c = ra.get().get();
		ra.bck( 0 );
		ra.bck( 0 );
		double d = ra.get().get();
		ra.fwd( 0 );
		ra.fwd( 1 );

		hessianMatrix[0][1] = hessianMatrix[1][0] =
				(
						(a - b) / 2
						-
						(c - d) / 2
				) / 2;
				/*(
						(img.get(x + 1, y + 1, z) - img.get(x - 1, y + 1, z)) / 2
						-
						(img.get(x + 1, y - 1, z) - img.get(x - 1, y - 1, z)) / 2
				) / 2;*/

		// xz
		ra.fwd( 0 );
		ra.fwd( 2 );
		a = ra.get().get();
		ra.bck( 0 );
		ra.bck( 0 );
		b = ra.get().get();
		ra.fwd( 0 );
		ra.fwd( 0 );
		ra.bck( 2 );
		ra.bck( 2 );
		c = ra.get().get();
		ra.bck( 0 );
		ra.bck( 0 );
		d = ra.get().get();
		ra.fwd( 0 );
		ra.fwd( 2 );

		hessianMatrix[0][2] = hessianMatrix[2][0] =
				(
						(a - b) / 2
						-
						(c - d) / 2
				) / 2;
				/*(
						(img.get(x + 1, y, z + 1) - img.get(x - 1, y, z + 1)) / 2
						-
						(img.get(x + 1, y, z - 1) - img.get(x - 1, y, z - 1)) / 2
				) / 2;*/

		// yz
		ra.fwd( 1 );
		ra.fwd( 2 );
		a = ra.get().get();
		ra.bck( 1 );
		ra.bck( 1 );
		b = ra.get().get();
		ra.fwd( 1 );
		ra.fwd( 1 );
		ra.bck( 2 );
		ra.bck( 2 );
		c = ra.get().get();
		ra.bck( 1 );
		ra.bck( 1 );
		d = ra.get().get();
		ra.fwd( 1 );
		ra.fwd( 2 );

		hessianMatrix[1][2] = hessianMatrix[2][1] =
				(
						(a - b) / 2
						-
						(c - d) / 2
				) / 2;
				/*(
						(img.get(x, y + 1, z + 1) - img.get(x, y - 1, z + 1)) / 2
						-
						(img.get(x, y + 1, z - 1) - img.get(x, y - 1, z - 1)) / 2
				) / 2;*/
	}
}
