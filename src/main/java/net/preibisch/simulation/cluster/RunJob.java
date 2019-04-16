package net.preibisch.simulation.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.preibisch.simulation.SimulateMultiViewAberrations;

public class RunJob
{
	public static boolean isCluster = false;

	public static void display( final RandomAccessibleInterval< FloatType > image, final String name )
	{
		if ( !isCluster )
		{
			ImagePlus imp = ImageJFunctions.wrapFloat( image, name ).duplicate();
			imp.setDimensions( 1, (int)image.dimension( 2 ), 1 );
			imp.show();
		}
	}

	public static void main( String[] args )
	{
		isCluster = true;

		final String dir = args[ 0 ].trim();
		final int zPlane = Integer.parseInt( args[ 1 ] );
		final double lsMiddle = Double.parseDouble( args[ 2 ] );
		final double lsEdge = Double.parseDouble( args[ 3 ] );
		final double ri = Double.parseDouble( args[ 4 ] );
		final boolean illum = Boolean.parseBoolean( args[ 5 ] ); // true > bottom, false > top

		//final double lsMiddle = 1.0;
		//final double lsEdge = 3.0;
		//final double ri = 1.005;

		final ExecutorService service = Executors.newFixedThreadPool( 1 );

		System.out.println( "dir=" + dir );
		System.out.println( "zPlane=" + zPlane );
		System.out.println( "lsMiddle=" + lsMiddle );
		System.out.println( "lsEdge=" + lsEdge );
		System.out.println( "ri=" + ri );
		System.out.println( "illum=" + illum );

		SimulateMultiViewAberrations.simulate( illum, lsMiddle, lsEdge, ri, dir, service, zPlane );
	}
}
