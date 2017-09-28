/*-
 * #%L
 * Library for simulating a multi-view acquisition including
 * attenuation, convolution, reduced sampling and poission noise.
 * %%
 * Copyright (C) 2014 - 2017 Multiview Simulation developers.
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
package net.preibisch.simulation.uncommons;
//============================================================================
//Copyright 2006-2010 Daniel W. Dyer
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//============================================================================

import java.util.Random;

/**
* Discrete random sequence that follows a
* <a href="http://en.wikipedia.org/wiki/Poisson_distribution" target="_top">Poisson
* distribution</a>.
* @author Daniel Dyer
*/
public class PoissonGenerator implements NumberGenerator<Integer>
{
private final Random rng;
private final NumberGenerator<Double> mean;


/**
* <p>Creates a generator of Poisson-distributed values.  The mean is
* determined by the provided {org.uncommons.maths.number.NumberGenerator}.  This means that
* the statistical parameters of this generator may change over time.
* One example of where this is useful is if the mean generator is attached
* to a GUI control that allows a user to tweak the parameters while a
* program is running.</p>
* <p>To create a Poisson generator with a constant mean, use the
* {@link #PoissonGenerator(double, Random)} constructor instead.</p>
* @param mean A {@link NumberGenerator} that provides the mean of the
* Poisson distribution used for the next generated value.
* @param rng The source of randomness.
*/
public PoissonGenerator(NumberGenerator<Double> mean,
                       Random rng)
{
   this.mean = mean;
   this.rng = rng;
}


/**
* Creates a generator of Poisson-distributed values from a distribution
* with the specified mean.
* @param mean The mean of the values generated.
* @param rng The source of randomness.
*/
public PoissonGenerator(double mean,
                       Random rng)
{
   this(new ConstantGenerator<Double>(mean), rng);
   if (mean <= 0)
   {
       throw new IllegalArgumentException("Mean must be a positive value.");
   }
}


/**
* {@inheritDoc}
*/
public Integer nextValue()
{
   int x = 0;
   double t = 0.0;
   while (true)
   {
       t -= Math.log(rng.nextDouble()) / mean.nextValue();
       if (t > 1.0)
       {
           break;
       }
       ++x;
   }
   return x;
}
}
