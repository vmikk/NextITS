#!/usr/bin/env -S awk -f

## Classify reads by the arrangement of primer sites (output of `sassy search`)

BEGIN {
  FS = "\t"
  OFS = "\t"
  print "SeqID", "Layout", "Sites", "Category"
  read = ""
}

function classify(    a, b) {
  if (n == 1) return "one_primer"
  if (n == 2) {
    a = lay[1]; b = lay[2]
    if (a == "F+" && b == "R-") return "valid_fwd"
    if (a == "R+" && b == "F-") return "valid_rev"
    if (a == "F+" && b == "F-") return "single_primer_F"
    if (a == "R+" && b == "R-") return "single_primer_R"
    if (substr(a, 1, 1) != substr(b, 1, 1)) return "misoriented"
    return "multiprimer"
  }
  return "multiprimer"
}

function flush(    i, layout, sites) {
  if (read == "") return
  layout = ""; sites = ""
  for (i = 1; i <= n; i++) {
    layout = layout lay[i]
    sites  = sites (i > 1 ? "," : "") st[i] "-" en[i]
  }
  print read, layout, sites, classify()
}

{
  if ($1 != read) {
    flush()
    read = $1
    n = 0
    split("", maxend)   # last site end per primer+strand
    split("", lastsite) # index of that site in lay/st/en
  }

  key = $2 $3
  start = $4 + 0
  end   = $5 + 0

  ## Overlaps an existing site of the same primer and strand -> extend that site
  if ((key in maxend) && start < maxend[key]) {
    if (end > maxend[key]) {
      maxend[key] = end
      en[lastsite[key]] = end
    }
    next
  }

  n++
  lay[n] = key
  st[n]  = start
  en[n]  = end
  maxend[key]   = end
  lastsite[key] = n
}

END { flush() }
