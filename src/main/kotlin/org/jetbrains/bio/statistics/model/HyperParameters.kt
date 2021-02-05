package org.jetbrains.bio.statistics.model

/**
 * A common interface for hyper parameters of Bayesian models.
 *
 * @author Sergei Lebedev
 * @since 02/03/15
 */
interface HyperParameters<T : HyperParameters<T, InstanceModel>,
        InstanceModel : ClassificationModel> {
    /**
     * Calculates Kullback-Leibler divergence between the current parameter
     * distribution and the given one.
     *
     * @see org.jetbrains.bio.statistics.distribution.KullbackLeibler
     */
    fun divergence(other: T): Double

    fun copy(): T

    /**
     * Returns log density of the given frequentist model parameters
     * under the current parameter distribution.
     */
    fun logDensity(model: InstanceModel): Double {
        throw UnsupportedOperationException()
    }

    fun sample(): InstanceModel = throw UnsupportedOperationException()
}
