package org.jetbrains.bio.statistics

import com.google.common.collect.Lists
import org.apache.commons.math3.exception.DimensionMismatchException
import org.apache.commons.math3.ml.clustering.CentroidCluster
import org.apache.commons.math3.ml.clustering.DoublePoint
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer
import org.apache.commons.math3.ml.clustering.MultiKMeansPlusPlusClusterer
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.dataframe.DataFrame
import java.util.*

/**
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 23/08/13
 */
object Clustering {
    /**
     * Clusters a [DataFrame] by rows using KMeans++.
     */
    fun kmeansPlusPlus(df: DataFrame, numClusters: Int, numTrials: Int = 1):
            List<CentroidCluster<DoublePoint>> {
        val points = ArrayList<DoublePoint>(df.rowsNumber)
        for (i in 0 until df.rowsNumber) {
            val row = df.rowAsDouble(i)
            if (row.asSequence().any(Double::isNaN)) {
                continue  // Silently skip NaNs.
            }

            points.add(DoublePoint(row))
        }

        val clusterer = KMeansPlusPlusClusterer<DoublePoint>(numClusters)
        return if (numTrials == 1) {
            clusterer.cluster(points)
        } else {
            MultiKMeansPlusPlusClusterer(clusterer, numTrials).cluster(points)
        }
    }

    /**
     * Computes summary statistics over each cluster from a given list.
     *
     * @param clusters a list of clusters to summarize.
     * @return within-cluster summaries, sorted by mean +- SE.
     */
    fun summarize(clusters: List<CentroidCluster<DoublePoint>>): List<SummaryStatistics> {
        val summaries = Lists.newArrayListWithCapacity<SummaryStatistics>(clusters.size)
        for (cluster in clusters) {
            val summary = SummaryStatistics()
            for (point in cluster.points) {
                val coordinates = point.point
                if (coordinates.size != 1) {
                    throw DimensionMismatchException(coordinates.size, 1)
                }

                summary.addValue(coordinates[0])
            }

            summaries.add(summary)
        }

        return summaries.sortedBy {it.mean + it.standardDeviation / FastMath.sqrt(it.n.toDouble())}
    }
}
