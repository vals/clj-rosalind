(ns clj-rosalind.core
    (:require [clojure.string :as string]))

(defn count-bases
  "Count the bases in a DNA string. Prints the counts in alphabetical order.
  http://rosalind.info/problems/dna/"
  [s]
  (string/join " " (vals (sort (frequencies s)))))

(defn transcribe-rna
  "Transcribes a DNA string to an RNA string.
  http://rosalind.info/problems/rna/"
  [s]
  (string/replace s #"T" "U"))

(def cb {"A" "T", "C" "G", "G" "C", "T" "A"})

(defn complement-dna
  "Return the complementary DNA string.
  http://rosalind.info/problems/revc/"
  [s]
  (string/replace s #"A|C|G|T" cb))

(defn reverse-complement
  "Returns reverse complement of a DNA string.
  http://rosalind.info/problems/revc/"
  [s]
  (string/reverse (complement-dna s)))

(defn fasta-to-map
  "Takes a string in FASTA format and converts it to a map of label keys
  and sequence values."
  [s]
  (into {} (remove #(< (count %) 2)
                   (map #(string/split %  #"\n" 2)
                        (string/split s #">")))))

(defn clean-sequence
  "Return a DNA string where anything not a base have been removed."
  [s]
  (apply str (filter #(#{\A,\C,\G,\T} %) s)))

(defn gc-content
  "Retrun the GC content of a DNA string."
  [s]
  (float (* 100 (/ (count (filter #(#{\C,\G} %) (clean-sequence s))) (count (clean-sequence s))))))

(defn update-values [m f & args]
  (reduce (fn [r [k v]] (assoc r k (apply f v args))) {} m))

(defn fasta-gc-content
  "Return the GC content of fasta records."
  [s]
  (update-values (fasta-to-map s) gc-content))

(defn top-fasta-gc-content
  "Return FASTA record with highest GC content

  http://rosalind.info/problems/gc/"
  [s]
  (apply max-key val (fasta-gc-content s)))
