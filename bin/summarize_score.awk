BEGIN {
    FS = "[[:space:]]+"
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "SCORE1_AVG") score_col = i
    }
    if (!score_col) {
        print "SCORE1_AVG column not found" > "/dev/stderr"
        exit 2
    }
    next
}

{
    value = $score_col + 0
    samples++
    sum += value
    sumsq += value * value
    if (samples == 1 || value < minimum) minimum = value
    if (samples == 1 || value > maximum) maximum = value
}

END {
    if (samples == 0) {
        print "No scored samples found" > "/dev/stderr"
        exit 2
    }
    mean = sum / samples
    variance = (sumsq / samples) - (mean * mean)
    if (variance < 0 && variance > -1e-30) variance = 0
    sd = sqrt(variance)
    missing = requested - matched
    status = (fraction < warn) ? "warning_low_variant_match" : "pass"
    print "trait\trequested_variants\tmatched_variants\tmissing_variants\tmatch_fraction\tsamples\tmean\tsd\tmin\tmax\tstatus"
    printf "%s\t%d\t%d\t%d\t%.8f\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%s\n", \
        trait, requested, matched, missing, fraction, samples, mean, sd, minimum, maximum, status
}
