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

/**
 * Interface for providing different types of sequences of numbers.  This is
 * a simple but powerful abstraction that provides considerable flexibility
 * in implementing classes that require numeric configuration.  Refer to the
 * implementations in this package for examples of how it can be used.
 * @param <T> The type (Integer, Long, Double, etc.) of number to generate.
 * @author Daniel Dyer
 */
public interface NumberGenerator<T extends Number>
{
    /**
     * @return The next value from the generator.
     */
    T nextValue();
}
