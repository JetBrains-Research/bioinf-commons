package org.jetbrains.bio.util

import org.junit.rules.TestRule
import org.junit.runner.Description
import org.junit.runners.model.Statement
import java.util.concurrent.atomic.AtomicInteger

/**
 * See https://wiki.saucelabs.com/display/DOCS/How+to+Deal+with+Flaky+Java+Tests
 */
class RetryRule(retries: Int) : TestRule {
    private val retryCount: AtomicInteger = AtomicInteger(retries)

    override fun apply(base: Statement, description: Description): Statement {
        return object : Statement() {
            override fun evaluate() {
                var caughtThrowable: Throwable? = null

                while (retryCount.getAndDecrement() > 0) {
                    try {
                        base.evaluate()
                        return
                    } catch (t: Throwable) {
                        if (description.getAnnotation(Retry::class.java) != null) {
                            if (retryCount.get() > 0) {
                                caughtThrowable = t
                                System.err.println("${description.displayName}: Failed, $retryCount retries remain")
                            } else {
                                throw caughtThrowable!!
                            }
                        } else {
                            // this test shouldn't be retried
                            throw t
                        }
                    }
                }
            }
        }
    }
}

@Retention(AnnotationRetention.RUNTIME)
annotation class Retry