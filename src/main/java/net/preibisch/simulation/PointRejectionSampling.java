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
