package org.jetbrains.bio.experiment

import org.jetbrains.bio.Configuration
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempFile
import org.junit.Test

class DataConfigExperimentTest {


    @Test(expected = IllegalStateException::class)
    fun differentConfig() {
        withTempFile("test1", ".bed") { p ->
            val k4me1 = DataConfig.load("""genome: to1
tracks:
   H3K4me1:
      t:
      - $p
""".reader(), "k4me1")
            val k4me3 = DataConfig.load("""genome: to1
tracks:
   H3K4me3:
      t:
      - $p
""".reader(), "k4me3")

            val dir = Configuration.experimentsPath
            (dir / "configs" / "k4me1" / "experiment" / "k4me1.yaml").bufferedWriter().use {
                k4me3.save(it)
            }
            // Check different config error
            object : DataConfigExperiment("experiment", k4me1) {
                override fun doCalculations() {
                    // Ignore
                }
            }.run()
        }
    }
}
