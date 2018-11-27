package net.preibisch.simulation.raytracing;

import java.util.ArrayList;
import java.util.Collection;

import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NoninvertibleModelException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;

/**
 * @author Stephan Preibisch and Varun Kapoor
 */
public class Lightsheet
{
	final int minNumPoints = 3;
	double a, b, c; // a*x*x + b*x + c

	public Lightsheet(
			final double center,
			final double thicknessCenter,
			final double length,
			final double thickNessEdges )
	{
		final ArrayList< Point > points = new ArrayList<>();
		points.add( new Point( new double[] { center, thicknessCenter } ) );
		points.add( new Point( new double[] { center - length/2, thickNessEdges } ) );
		points.add( new Point( new double[] { center + length/2, thickNessEdges } ) );

		try
		{
			fit( points );

			System.out.println( center - length/2 + ": " + predict( center - length/2 ) );
			System.out.println( center + ": " + predict( center ) );
			System.out.println( center + length/2 + ": " + predict( center + length/2 ) );
		}
		catch ( NotEnoughDataPointsException | IllDefinedDataPointsException e )
		{
			throw new RuntimeException( "Could not fit quadratic function for lightsheet: " + e );
		}
	}

	public Lightsheet( final double a, final double b, final double c )
	{
		this.a = a;
		this.b = b;
		this.c = c;
	}

	public double getA(){ return a; }
	public double getB(){ return b; }
	public double getC(){ return c; }

	public double predict( final double x ) { return a*x*x + b*x + c; }

	public void fit( final Collection< Point > points ) throws NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		final int numPoints = points.size();

		if ( numPoints < minNumPoints )
			throw new NotEnoughDataPointsException( "Not enough points, at least " + minNumPoints + " are necessary and available are: " + numPoints );

		// compute matrices
		final double[] delta = new double[ 9 ];
		final double[] tetha = new double[ 3 ];
		
		for ( final Point p : points )
		{
			final double x = p.getW()[ 0 ];
			final double y = p.getW()[ 1 ];

			final double xx = x*x;
			final double xxx = xx*x;

			delta[ 0 ] += xx * xx;
			delta[ 1 ] += xxx;
			delta[ 2 ] += xx;

			delta[ 3 ] += xxx;
			delta[ 4 ] += xx;
			delta[ 5 ] += x;

			delta[ 6 ] += xx;
			delta[ 7 ] += x;
			delta[ 8 ] += 1;

			tetha[ 0 ] += xx * y;
			tetha[ 1 ] += x * y;
			tetha[ 2 ] += y;
		}

		// invert matrix
		try
		{
			invert3x3( delta );
		}
		catch ( final NoninvertibleModelException e )
		{
			this.a = this.b = this.c = 0;
			throw new IllDefinedDataPointsException( "Cannot not invert Delta-Matrix, failed to fit function" );
		}

		this.a = delta[ 0 ] * tetha[ 0 ] + delta[ 1 ] * tetha[ 1 ] + delta[ 2 ] * tetha[ 2 ];
		this.b = delta[ 3 ] * tetha[ 0 ] + delta[ 4 ] * tetha[ 1 ] + delta[ 5 ] * tetha[ 2 ];
		this.c = delta[ 6 ] * tetha[ 0 ] + delta[ 7 ] * tetha[ 1 ] + delta[ 8 ] * tetha[ 2 ];
	}

	final static public double det3x3( final double[] a )
	{
		assert a.length == 9 : "Matrix3x3 supports 3x3 double[][] only.";
		
		return
			a[ 0 ] * a[ 4 ] * a[ 8 ] +
			a[ 3 ] * a[ 7 ] * a[ 2 ] +
			a[ 6 ] * a[ 1 ] * a[ 5 ] -
			a[ 2 ] * a[ 4 ] * a[ 6 ] -
			a[ 5 ] * a[ 7 ] * a[ 0 ] -
			a[ 8 ] * a[ 1 ] * a[ 3 ];
	}

	final static public void invert3x3( final double[] m ) throws NoninvertibleModelException
	{
		assert m.length == 9 : "Matrix3x3 supports 3x3 double[][] only.";
		
		final double det = det3x3( m );
		if ( det == 0 ) throw new NoninvertibleModelException( "Matrix not invertible." );
		
		final double i00 = ( m[ 4 ] * m[ 8 ] - m[ 5 ] * m[ 7 ] ) / det;
		final double i01 = ( m[ 2 ] * m[ 7 ] - m[ 1 ] * m[ 8 ] ) / det;
		final double i02 = ( m[ 1 ] * m[ 5 ] - m[ 2 ] * m[ 4 ] ) / det;
		
		final double i10 = ( m[ 5 ] * m[ 6 ] - m[ 3 ] * m[ 8 ] ) / det;
		final double i11 = ( m[ 0 ] * m[ 8 ] - m[ 2 ] * m[ 6 ] ) / det;
		final double i12 = ( m[ 2 ] * m[ 3 ] - m[ 0 ] * m[ 5 ] ) / det;
		
		final double i20 = ( m[ 3 ] * m[ 7 ] - m[ 4 ] * m[ 6 ] ) / det;
		final double i21 = ( m[ 1 ] * m[ 6 ] - m[ 0 ] * m[ 7 ] ) / det;
		final double i22 = ( m[ 0 ] * m[ 4 ] - m[ 1 ] * m[ 3 ] ) / det;
		
		m[ 0 ] = i00;
		m[ 1 ] = i01;
		m[ 2 ] = i02;

		m[ 3 ] = i10;
		m[ 4 ] = i11;
		m[ 5 ] = i12;

		m[ 6 ] = i20;
		m[ 7 ] = i21;
		m[ 8 ] = i22;
	}

}
