(ns clj-rosalind.core)

(defn count-bases
  "Count the bases in a DNA string. Prints the counts in alphabetical order.
  http://rosalind.info/problems/dna/"
  [s]
  (clojure.string/join " " (vals (sort (frequencies s)))))

(defn transcribe-rna
    "Transcribes a DNA string to an RNA string.
    http://rosalind.info/problems/rna/"
    [s]
    (clojure.string/replace s #"T" "U"))

(def cb {"A" "T", "C" "G", "G" "C", "T" "A"})

(defn complement-dna
    "Return the complementary DNA string.
    http://rosalind.info/problems/revc/"
    [s]
    (clojure.string/replace s #"A|C|G|T" cb))

(defn reverse-complement
    "Returns reverse complement of a DNA string.
    http://rosalind.info/problems/revc/"
    [s]
    (clojure.string/reverse (complement-dna s)))
