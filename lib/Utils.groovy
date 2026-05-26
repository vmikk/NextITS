/*
  Utility functions for NextITS pipeline.
  Placed in lib/ so they are available on the JVM classpath before config files are evaluated,
  which is required by the Nextflow v26 strict syntax parser
 (function definitions are not allowed in .config files)
*/
class Utils {

    // Ensure that a resource requirement does not exceed a configured maximum
    // Call from config closures as: { Utils.checkMax(3.h * task.attempt, 'time', params) }
    static def checkMax(obj, type, params) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min(obj, params.max_cpus as int)
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }

    // Scale memory requirements based on retry attempt number.
    //   1st attempt = 1.0x
    //   2nd         = 1.25x
    //   3rd         = 1.5x
    //   4th         = 1.75x
    //   max cap     = 3.0x
    // Call from config closures as: { Utils.scaleRetryMemory(8.GB, task.attempt, params) }
    static def scaleRetryMemory(base, attempt, params) {
        def retries = Math.max(attempt as int, 1) - 1
        def factor  = Math.min(1.0 + (retries * 0.25), 3.0)
        def scaledBytes = (base.toBytes() * factor) as long
        return checkMax(scaledBytes.B, 'memory', params)
    }

}
