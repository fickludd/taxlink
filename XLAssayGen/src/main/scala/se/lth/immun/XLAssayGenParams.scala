package se.lth.immun

import java.io.File
import se.jt.Params
import se.lth.immun.xlink._
import se.lth.immun.chem.EPeptideFragment

class XLAssayGenParams extends Params {

	import Params._
	
	val xlMaster = ReqString("Input xl master format file")
	val spectralLib = ""			## "Path to spectral library to use for fragment selection (NOT IMPLEMENTED YET)"
	val noCarbAmidoMethyl = false	## "Set to disable fixed carbamidomethylation of cysteines in Kojak XLs"
	val isotopeOccurenceMass = 0.9	## "Take top occuring precursor isotopes until their summed occurence is above this"
	val qValue = 1.0				## "Set maximum q-value allowed for rows with perculator result"
	val useIncompleteRows = false	## "Use rows that lack taxlink or kojak/percolator results"
	val xlDeuteriums		= 12 	## "The number of deuteriums in the labeled xlinkers"
	val fragments = "y,b"			## "Fragment types to generate [comma-sep list of a,b,c,x,y,z]"
	
	val outHumanImg = true 			## "Set to enable humanreadable output figures"
	val outDir			= ""		## "output directory (by default same as input feature file)"
	val outName			= ""		## "basename for output files (by default same as input feature file)"
	
	def outBase = {
		val taxFile = new File(xlMaster)
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
		stripExt(stripExt(path, ".tsv"), ".csv")
		
	lazy val kojakConf = 
		if (noCarbAmidoMethyl) KojakXLink.Conf(XLink.DSS, Map())
		else KojakXLink.DSS_CARB_CONF
		
	var FRAG_TYPES:Array[EPeptideFragment] = _
	def parseFragTypes = 
		FRAG_TYPES = fragments.value.split(",").map(str => EPeptideFragment.fromChar(str.head)).toArray
	
}