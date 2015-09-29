package se.lth.immun

import java.io.File
import se.jt.Params

class TaxPerKojakParams extends Params {

	import Params._
	
	val taxlink = ""	## "Taxlink output file"
	val kojak = ""		## "Kojak output file"
	val perc = ""		## "Percolator (post kojak) output file"
	val mzML = ""		## "mzML / mzML.gz file used for parsing retention time"
	
	val qValue = 1.0	## "FDR threshold to apply on kojak q-values"
	val matchPPM = 3.0	## "PPM threshold for taxlink feature - kojak matching"
	val matchRT = 0.5	## "Allowed rt diff between kojak ID and taxlink average apex"
	val excludeProteinsWith = ""	## "Exclude any XLs for proteins containing this tag in their protein ID"
	val outDir			= ""		## "output directory (by default same as input feature file)"
	val outName			= ""		## "basename for output files (by default same as input feature file)"
	val verbose = false			## "Set to get more output"
	
	def outBase = {
		val taxFile = new File(taxlink)
		val dir = 
			if (outDir.value != "") outDir.value
			else taxFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExts(taxFile.getName)
		(dir, name) 
	}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(stripExt(path, ".tsv"), ".csv"), ".isopairs")
	
	def taxOpt = 
		if (taxlink.value == "") None
		else Some(taxlink.value)
	
	def percOpt = 
		if (perc.value == "") None
		else Some(perc.value)
		
	def kojakOpt = 
		if (kojak.value == "") None
		else Some(kojak.value)
		
	def mzMLOpt = 
		if (mzML.value == "") None
		else Some(mzML.value)
	
	val outSep = "\t"
	val outQuote = "\""
}