package net.preibisch.simulation.raytracing;

public class ClearingMap implements RefractiveIndexMap
{

	@Override
	public double getRI( final double intensity )
	{
		return 0;
	}

}
