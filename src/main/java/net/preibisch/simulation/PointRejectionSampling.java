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
package net.preibisch.simulation;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.type.numeric.RealType;

public class PointRejectionSampling
{
	public static <T extends RealType< T > > List<RealPoint> sampleRealPoints(RealInterval interval, int nSamples, RealRandomAccessible< T > density, Random rnd)
	{
		final List<RealPoint> samples = new ArrayList<>();
		final int nDim = interval.numDimensions();
		final RealRandomAccess< T > ra = density.realRandomAccess();
		while(samples.size() < nSamples)
		{
			double[] pos = new double[nDim];
			for (int d=0; d<nDim; d++)
			{
				pos[d] = interval.realMin( d ) + rnd.nextDouble() * (interval.realMax( d ) - interval.realMin( d ));
			}
			double p = rnd.nextDouble();
			ra.setPosition( pos );
			if (p < ra.get().getRealDouble())
				samples.add( new RealPoint( pos ) );
		}
		return samples;
	}
	
}
