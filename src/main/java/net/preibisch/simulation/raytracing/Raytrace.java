/*-
 * #%L
 * Library for simulating a multi-view acquisition including
 * attenuation, convolution, reduced sampling and poission noise.
 * %%
 * Copyright (C) 2014 - 2026 Multiview Simulation developers.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
package net.preibisch.simulation.raytracing;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Raytrace
{
	public static void reflect( final double[] i, final double[] n, final double[] r )
	{
		// r = i−2(i·n)n
		final double dotP = i[0]*n[0] + i[1]*n[1] + i[2]*n[2];

		r[ 0 ] = i[ 0 ] - 2*dotP*n[ 0 ];
		r[ 1 ] = i[ 1 ] - 2*dotP*n[ 1 ];
		r[ 2 ] = i[ 2 ] - 2*dotP*n[ 2 ];
	}

	public static double refract( final double[] i, final double[] n, final double n0, final double n1, final double thetaI, final double[] t )
	{
		final double deltaN = n0 / n1;
		final double thetaT = Math.asin( deltaN * Math.sin( thetaI ) );

		// total internal reflection
		if ( Double.isNaN( thetaT ) )
			return thetaT;

		final double cosThetaI = Math.cos( thetaI );
		final double sinThetaT = Math.sin( thetaT );

		t[ 0 ] = deltaN * i[ 0 ] - n[ 0 ] * ( deltaN * cosThetaI - Math.sqrt( 1 - sinThetaT*sinThetaT ));
		t[ 1 ] = deltaN * i[ 1 ] - n[ 1 ] * ( deltaN * cosThetaI - Math.sqrt( 1 - sinThetaT*sinThetaT ));
		t[ 2 ] = deltaN * i[ 2 ] - n[ 2 ] * ( deltaN * cosThetaI - Math.sqrt( 1 - sinThetaT*sinThetaT ));

		return thetaT;
	}

	public static double length( final double[] v )
	{
		return Math.sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
	}

	public static void norm( final double[] v )
	{
		final double l = length( v );

		v[ 0 ] /= l;
		v[ 1 ] /= l;
		v[ 2 ] /= l;
	}

	public static double incidentAngle( final double[] i, final double[] n )
	{
		double thetaI = Math.acos( ( n[0]*i[0] + n[1]*i[1] + n[2]*i[2]) / ( Math.sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] ) * Math.sqrt( i[0]*i[0] + i[1]*i[1] + i[2]*i[2] ) ) );

		if ( thetaI >= Math.PI / 2 )
		{
			// invert normal vector & eigenvalue
			n[ 0 ] *= -1;
			n[ 1 ] *= -1;
			n[ 2 ] *= -1;
			//ev *= -1;

			// adjust angle
			thetaI -= Math.PI / 2;
			//dotP = Math.acos( ( nx*bx + ny*by + nz*bz) / ( Math.sqrt( nx*nx + ny*ny + nz*nz ) * Math.sqrt( bx*bx + by*by + bz*bz ) ) );
		}

		return thetaI;
	}


	public static void drawSimpleImage( final RandomAccessibleInterval< FloatType > randomAccessible, final float low, final float high )
	{
		final Cursor< FloatType > c = Views.iterable( randomAccessible ).localizingCursor();

		while( c.hasNext() )
		{
			c.fwd();
			boolean all = true;
			//for ( int d = 0; d < n; ++d )
			//	if ( c.getIntPosition( d ) > 100 || c.getIntPosition( d ) < 10 )
			//		all = false;

			if ( c.getIntPosition( 0 ) < randomAccessible.dimension( 0 ) / 2 )
			{
				if ( c.getIntPosition( 0 ) > c.getIntPosition( 1 ) )
					all = true;
				else
					all = false;
			}
			else
			{
				if ( randomAccessible.dimension( 0 ) - c.getIntPosition( 0 ) > c.getIntPosition( 1 ) )
					all = true;
				else
					all = false;
			}

			if ( all )
				c.get().set( high );
			else
				c.get().set( low );
		}
	}

}
