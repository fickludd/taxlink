package se.lth.immun

import se.jt.Params
import se.lth.immun.chem.Element
import se.lth.immun.xlink._
import java.io.File

class TaxlinkParams(val name:String, val version:String) extends Params {
	import Params._
	
	// USER PARAMS
	val featMzTolerance		= 0.1 	## "the maximum offset compared to theoretical mass that should be considered (Da)"
	val allowMissedMonoIsotope = false //## "Check if features have missed the smallest isotope"
	val xlDeuteriums		= 12 	## "The number of deuteriums in the labeled xlinkers"
	val xlRtDiff			= 5.0	## "the expected avg. rt offset of the light and heavy XLs (s)"
	val xlRtDiff95			= 5.0	## "the expected spread in rt offset of the light and heavy XLs (s)"
	val xlPairInt95 = (math.log(2) - math.log(1)) ## "the expected spread in logged intensity of the light and heavy XLs"
	val mzAccuracy95 		= 0.01 ## "the expected spread in measured vs. theoretical m/z (Da)"
	val targetRtDiff		= 1.0	## "Allowed difference of feature and target rt"
	val noCarbAmidoMethyl = false	## "Set to disable fixed carbamidomethylation of cysteines in Kojak XLs"
	//val gradientLength			= 8000.0 ## "the gradient length (s)"
	val outDir			= ""		## "output directory (by default same as input feature file)"
	val outName			= ""		## "basename for output files (by default same as input feature file)"
	
	val featurePath = ReqString("The dinosaur feature file to use")
	val xlPath = ReqString("File listing the sought Xlinks")
	
	
	lazy val kojakConf = 
		if (noCarbAmidoMethyl) KojakXLink.Conf(XLink.DSS, Map())
		else KojakXLink.DSS_CARB_CONF
	
	
	def xlinkerDeltaMass = 
		xlDeuteriums * (Element.H2.monoisotopicWeight - Element.H.monoisotopicWeight)
		
	def outBase = {
		val featFile = new File(featurePath)
		val dir = 
			if (outDir.value != "") outDir.value
			else featFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExts(featFile.getName)
		(dir, name) 
	}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(stripExt(path, ".tsv"), ".csv"), ".features")
		
	val outSep = "\t"
	val outQuote = "\""
}