package net.preibisch.simulation.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class CreateScripts
{
	public static void main( String[] args )
	{
		createScripts();
	}

	public static void createScripts()
	{
		PrintWriter f = openFileWrite( new File( "submitAll.sh" ) );

		for ( int z = 0; z < 289; ++z )
		{
			for ( int i = 0; i <= 1; ++i )
			{
				String illum = i == 0 ? "false" : "true";
				String jobFilename =  "./job_"+ z + "_" + illum + ".sh";
				String command = "./java -Xms4000M -Xmx4000M -jar multiview-simulation-0.2.1-SNAPSHOT-jar-with-dependencies.jar /fast/AG_Preibisch/Stephan/ " + z + " 3.0 1.0 1.1 " + illum;

				PrintWriter file = openFileWrite( new File( jobFilename ) );
				file.println( "#!/bin/sh" );
				file.println( command );
				file.close();

				f.println( "qsub -l h_vmem=34G " + jobFilename );
			}
		}

		f.close();
	}

	public static BufferedReader openFileRead(final File file)
	{
		BufferedReader inputFile;
		try
		{
			inputFile = new BufferedReader(new FileReader(file));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileRead(): " + e);
			inputFile = null;
		}
		return (inputFile);
	}

	public static PrintWriter openFileWrite(final File file)
	{
		PrintWriter outputFile;
		try
		{
			outputFile = new PrintWriter(new FileWriter(file));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileWrite(): " + e);
			outputFile = null;
		}
		return (outputFile);
	}

}
