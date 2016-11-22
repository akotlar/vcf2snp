package main

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	// "sort"
	"regexp"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
	"runtime/pprof"
)

func main() {
	inPath := flag.String("inPath", "", "The input .vcf file path")
	outPath := flag.String("outPath", "", "The output path (optional, default is stdout)")

	cpuprofile := flag.String("cpuProfile", "", "write cpu profile to file")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	iupac := map[string]string{
		"AA": "A", "AC": "M", "AG": "R", "AT": "W",
		"CA": "M", "CC": "C", "CG": "S", "CT": "Y",
		"GA": "R", "GC": "S", "GG": "G", "GT": "K",
		"TA": "W", "TC": "Y", "TG": "K", "TT": "T",
		"A": "A", "C": "C", "G": "G", "T": "T",
		"D": "D", "I": "I", "DD": "D", "II": "I",
	}

	inFh := (*os.File)(nil)

	// make sure it gets closed
	defer inFh.Close()

	if *inPath != "" {
		var err error
		inFh, err = os.Open(*inPath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	outFh := (*os.File)(nil)

	if *outPath != "" {
		var err error
		outFh, err = os.OpenFile(*outPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		outFh = os.Stdout
	}

	defer outFh.Close()

	rdr, err := vcfgo.NewReader(inFh, false)

	if err != nil {
		panic(err)
	}

	var out bytes.Buffer

	numSamples := len(rdr.Header.SampleNames)
	out.WriteString("Fragment\tPosition\tReference\tType\tAlleles\tAllele_Counts\t")
	for index, name := range rdr.Header.SampleNames {
		out.WriteString(name)

		if index-1 < numSamples {
			out.WriteString("\t")
		}
	}

	out.WriteString("\n")

	outFh.WriteString(out.String())

	// TODO: VCF files may put multiple alleles on multiple lines, same position
	// I vote that we just leave that alone
	// Our output reflects your input: GIGO
	// Else we'll need to gather those lines, slower, greater possibility of errors
	insCount := 0
	delCount := 0

	var chrBuffer bytes.Buffer
	var gt bytes.Buffer
	var sampleOut bytes.Buffer
	var altBuff bytes.Buffer
	var altCountBuff bytes.Buffer
	var gtIndex int

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}

		// TODO: modify the source library to allow auto appending chr
		if !strings.Contains(variant.Chromosome, "chr") {
			// var chrBuffer bytes.Buffer
			chrBuffer.Reset()

			chrBuffer.WriteString("chr")
			chrBuffer.WriteString(variant.Chromosome)
			variant.Chromosome = chrBuffer.String()
		}

		pos, ref, alt, alleleTypes, err := getCall(variant.Pos, variant.Ref(), variant.Alt())

		if err != nil {
			log.Print(err)
			continue
		}

		var callType string

		// fmt.Printf("Original: %s\t%d\t%s\t%s\n", variant.Chromosome, variant.Pos, variant.Ref(), variant.Alt())
		// fmt.Printf("New: %s\t%d\t%s\t%s\t%s\n", variant.Chromosome, pos, ref, alt, callType)

		// var out bytes.Buffer
		out.Reset()

		out.WriteString(variant.Chromosome)
		out.WriteString("\t")
		out.WriteString(strconv.FormatUint(pos, 10))
		out.WriteString("\t")
		out.WriteString(ref)

		// var sampleOut bytes.Buffer
		sampleOut.Reset()

		/************** Accumulate our SNP format genotypes: **********************
												genotype\tconfidence
												EX: R\t1
												For now treat all genotype confidence as 1
												Since VCF doesn't require a GQ call per sample
		***************************************************************************/
		seenGenos := make(map[string]int)

		for _, sample := range variant.Samples {
			delCount = 0
			insCount = 0

			// var gt bytes.Buffer
			gt.Reset()

			if len(sample.GT) > 2 {
				//TODO: Just skip the line?
				log.Fatal("Sample more than diploid, input file is weird")
			}

			for _, genotype := range sample.GT {
				if genotype == 0 {
					gt.WriteString(ref)
					if _, exists := seenGenos[ref]; !exists {
						seenGenos[ref] = 1
					} else {
						seenGenos[ref]++
					}

				} else {
					gtIndex = genotype - 1

					if alleleTypes[gtIndex] == "D" {
						delCount++
					} else if alleleTypes[gtIndex] == "I" {
						insCount++
					} else {
						gt.WriteString(alt[gtIndex])
					}

					if _, exists := seenGenos[alt[gtIndex]]; !exists {
						seenGenos[alt[gtIndex]] = 1
					} else {
						seenGenos[alt[gtIndex]]++
					}
				}
			}

			if delCount > 0 {
				if delCount == len(sample.GT) {
					sampleOut.WriteString("D")
				} else {
					sampleOut.WriteString("E")
				}
			} else if insCount > 0 {
				if insCount == len(sample.GT) {
					sampleOut.WriteString("I")
				} else {
					sampleOut.WriteString("H")
				}
			} else {
				// Note that if we left out delCount and insCount we could only get
				// homozygous indels
				sampleOut.WriteString(iupac[gt.String()])
			}

			// Fake quality
			sampleOut.WriteString("\t1\t")
		}

		altBuff.Reset()
		altCountBuff.Reset()

		seenRef := 0
		genoCount := 0
		for geno, count := range seenGenos {
			if seenRef == 0 && geno == ref {
				seenRef = 1
			}

			altBuff.WriteString(geno)
			altCountBuff.WriteString(strconv.Itoa(count))

			if genoCount++; genoCount != len(seenGenos) {
				altBuff.WriteString(",")
				altCountBuff.WriteString(",")
			}
		}

		// if len(seenGenos) == seenRef {
		// 	log.Printf("Saw only 1 genotype, which was ref (%s), skipping row %s:%s",
		// 		ref, variant.Chromosome, strconv.FormatUint(pos, 10))
		// 	continue
		// }

		if len(seenGenos)-seenRef > 1 {
			callType = "MULTIALLELIC"
		} else if delCount > 0 {
			callType = "DEL"
		} else if insCount > 0 {
			callType = "INS"
		} else {
			callType = "SNP"
		}

		out.WriteString("\t")
		out.WriteString(callType)
		out.WriteString("\t")
		out.WriteString(altBuff.String())
		out.WriteString("\t")
		out.WriteString(altCountBuff.String())
		out.WriteString("\t")
		out.WriteString(strings.TrimRight(sampleOut.String(), "\t"))
		out.WriteString("\n")

		outFh.WriteString(out.String())
	}

	log.Print(rdr.Error())

	pprof.StopCPUProfile()
}

var g = regexp.MustCompile(`^[ACTG]+$`)

func getCall(position uint64, reference string, alternateAlleles []string) (uint64, string, []string, []string, error) {
	var pos uint64

	// Right normalize the reference
	// If deletion, VCF reports the allele as 1 in length
	// and the ref as N length, with the length of the deletion being N - 1
	// SNP is the opposite on both counts, it's right normalized, meaning
	// the Ref is always 1 base, and the Alt is -N where N is the number of
	// bases deleted, inclusive of the reference
	ref := string(reference[len(reference)-1])
	alt := make([]string, 0, len(alternateAlleles))
	alleleTypes := make([]string, 0, 3)

	if !g.MatchString(reference) {
		return pos, ref, alt, alleleTypes, fmt.Errorf("Ref alleles must be composed of ACTG, got %s", reference)
	}

	for _, alternate := range alternateAlleles {
		if !g.MatchString(alternate) {
			return pos, ref, alt, alleleTypes, fmt.Errorf("Alt alleles must be composed of ACTG, got %s", alternate)
		}

		if len(reference) > len(alternate) {
			alleleTypes = append(alleleTypes, "D")

			sizeDiff := len(reference) - len(alternate)
			pos = position + uint64(sizeDiff)

			alt = append(alt, strconv.Itoa(-sizeDiff))
		} else if len(reference) < len(alternate) {
			alleleTypes = append(alleleTypes, "I")

			var buffer bytes.Buffer

			buffer.WriteString("+")
			// Seqant uses the string EXCLUDING the reference
			// as the allele, since it really makes 0 sense to include the reference
			// in the alternate string, under the definition alt == !ref
			buffer.WriteString(alternate[len(reference):])

			if len(buffer.String()) == 1 {
				return pos, ref, alt, alleleTypes, fmt.Errorf(`Failed to make alt for
					INS: %d, %s`, position, reference)
			}

			// Only deletions require us to modify our position
			if pos == 0 {
				pos = position
			}

			alt = append(alt, buffer.String())
		} else {
			alleleTypes = append(alleleTypes, "S")

			if pos == 0 {
				pos = position
			}

			alt = append(alt, alternate[len(reference)-1:])
		}
	}

	return pos, ref, alt, alleleTypes, nil
}
