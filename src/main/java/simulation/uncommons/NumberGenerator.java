package simulation.uncommons;

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
