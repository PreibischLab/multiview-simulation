package net.preibisch.simulation.cluster;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class RunJob
{
	public static boolean isCluster = false;

	public static void display( final RandomAccessibleInterval< FloatType > image, final String name )
	{
		if ( !isCluster )
			ImageJFunctions.wrapFloat( image, name ).show();;
	}
}
