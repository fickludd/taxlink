package se.lth.immun

import se.jt.CLIApp
import se.jt.CLIBar

import java.io.File
import java.io.FileReader
import java.io.Reader
import java.io.BufferedReader
import java.io.FileWriter
import java.io.BufferedWriter
import java.util.Properties

import scala.io.Source
import collection.mutable.ArrayBuffer
import collection.JavaConversions._
import collection.mutable.ArrayBuilder

import se.lth.immun.traml.ghost._
import se.lth.immun.xml._
import se.lth.immun.chem._
import se.lth.immun.xlink.XLink
import se.lth.immun.unimod.UniMod

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation

object Taxlink extends CLIApp {
	
	val ISOTOPE_PATTERN_DIFF = 1.00286864
		
	val nullFeature = new Feature(0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
	val nullScores = new Scores(0,0,0,0,0,0,0)
	
	
	case class XL(mass:Double, xlink:XLink, str:String, rt:Option[Double], z:Option[Int])
	case class FeaturePair(index1:Int, light:DinoFeature, heavy:DinoFeature)
	case class XLinkPairs(xl:XL, pairs:Seq[FeaturePair])
	
	case class DinoFeature(
		val mz:Double,
		val mostAbundant:Double,
		val z:Int,
		val rtStart:Double,
		val rtApex:Double,
		val rtEnd:Double,
		val nIsotopes:Int,
		val nScans:Int,
		val averagineCorr:Double,
		val mass:Double,
		val massCalib:Double,
		val intensityApex:Double,
		val intensitySum:Double
	)
	
	
	val cliBar = new CLIBar(30, true)
	val pc = new PearsonsCorrelation()
	
	var properties = new Properties
    properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.artifactId")
	val version 	= properties.getProperty("pom.version")

	implicit val params = new TaxlinkParams(name, version)	
	
	
	def main(args:Array[String]):Unit = {
		
		
		val t0 = System.currentTimeMillis
		
		failOnError(parseArgs(name, version, args, params, List("featurePath", "xlPath"), None))
    			
		val (outDir, outName) = params.outBase
		
		println("   "+name+" "+version)
		println("   dino feature file: " + params.featurePath.value)
    	println("          xlink file: " + params.xlPath.value)
    	println("             out dir: " + outDir)
    	println("            out name: " + outName)
		println
		println("reading features...")
		val features = readFeatures(new File(params.featurePath))
		
		println("reading xlinks...")
		val xlinks = readXlinks(new File(params.xlPath), params)
		
		println("matching xlinks to features...")
		val matches = matchFeatures(features, xlinks, params)
		
		println("writing %d matches".format(matches.map(_.pairs.length).sum))
		val outFile = new File(outDir, outName + ".isopairs.csv")
		writeMatches(matches, outFile)
		/*
		println("reading targets...")
		val targets = getTargets(gTraML)
		
		println("analyzing file...")
		val scoredFeatures = extractAndScore(features, targets)
		
		println("writing output file...")
		writeScoredFeatures(scoredFeatures)
		println("done")
		*/
		
		val t3 	= System.currentTimeMillis
		println("total time: "+niceTiming(t3-t0))
	}
	
	def readFeatures(f:File):ArrayBuffer[DinoFeature] = {
		var headerParsed = false
		
		var iMZ = -1
		var iMOST_ABUNDANT_MZ = -1
		var iZ = -1
		var iRT_START = -1
		var iRT_APEX = -1
		var iRT_END = -1
		var iN_ISO = -1
		var iN_SCAN = -1
		var iAVE_CORR = -1
		var iMASS = -1
		var iMASS_CAL = -1
		var iINT_APEX = -1
		var iINT_SUM = -1
		
		var res = new ArrayBuffer[DinoFeature]
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
				iMZ = cols.indexOf("mz")
				iMOST_ABUNDANT_MZ = cols.indexOf("mostAbundantMz")
				iZ = cols.indexOf("charge")
				iRT_START = cols.indexOf("rtStart")
				iRT_APEX = cols.indexOf("rtApex")
				iRT_END = cols.indexOf("rtEnd")
				iN_ISO = cols.indexOf("nIsotopes")
				iN_SCAN = cols.indexOf("nScans")
				iAVE_CORR = cols.indexOf("averagineCorr")
				iMASS = cols.indexOf("mass")
				iMASS_CAL = cols.indexOf("massCalib")
				iINT_APEX = cols.indexOf("intensityApex")
				iINT_SUM = cols.indexOf("intensitySum")
				headerParsed = true
			} else {
				res += DinoFeature(
					cols(iMZ).toDouble,
					cols(iMOST_ABUNDANT_MZ).toDouble,
					cols(iZ).toInt,
					cols(iRT_START).toDouble * 60,
					cols(iRT_APEX).toDouble * 60,
					cols(iRT_END).toDouble * 60,
					cols(iN_ISO).toInt,
					cols(iN_SCAN).toInt,
					cols(iAVE_CORR).toDouble,
					cols(iMASS).toDouble,
					cols(iMASS_CAL).toDouble,
					cols(iINT_APEX).toDouble,
					cols(iINT_SUM).toDouble
				)
			}
		}
		res
	}
	
	
	def readXlinks(f:File, params:TaxlinkParams):ArrayBuffer[XL] = {
		var headerParsed = false
		
		var iSEQ = -1
		var iRT = -1
		var iZ = -1
		
		var res = new ArrayBuffer[XL]
		for (line <- Source.fromFile(f).getLines) {
			val cols = line.split("\t")map(_.trim)
			if (!headerParsed) {
				iSEQ = cols.indexOf("sequence")
				iRT = cols.indexOf("rt")
				iZ = cols.indexOf("z")
				headerParsed = true
			} else {
				val str = cols(iSEQ)
				val rt = if (iRT < 0) None else Some(cols(iRT).toDouble)
				val z = if (iZ < 0) None else Some(cols(iZ).toInt)
				XLink.fromString(str, UniMod.parseUniModSequence, params.kojakConf) match {
					case Left(xlinkMol) =>
						res += XL(xlinkMol.monoisotopicMass, xlinkMol, str, rt, z)
					case Right(err) =>
						println(err)
				}
			}
		}
		res
	}
	
	
	
	def targetMatch(xl:XL, feature:DinoFeature) = 
		xl.rt.map(rt => 
			math.abs(rt - feature.rtApex) < params.targetRtDiff
		).getOrElse(true) &&
		xl.z.map(_ == feature.z).getOrElse(true)
		
	
	
	
	def matchFeatures(
			features:Seq[DinoFeature], 
			xlinks:Seq[XL], 
			params:TaxlinkParams
	):Seq[XLinkPairs] = {
		val fSorted = features.sortBy(_.mass)
		val xSorted = xlinks.flatMap(xl => 
						if (params.allowMissedMonoIsotope)
							Array(
								xl,
								XL(xl.mass + ISOTOPE_PATTERN_DIFF, xl.xlink, xl.str, xl.rt, xl.z)
							)
						else Array(xl)
					).sortBy(_.mass)
		
		var iStart = features.indexWhere(f => f.mz >= xSorted.head.mass - params.mzAccuracy95)
		if (iStart < 0) iStart = 0
		var iEnd = iStart
		for (xl <- xSorted) yield {
			while (iStart < fSorted.length && 
					fSorted(iStart).mass < xl.mass - params.featMzTolerance)
				iStart += 1
			
			iEnd = math.max(iStart, iEnd)
			while (iEnd < fSorted.length && 
					fSorted(iEnd).mass < xl.mass + params.featMzTolerance)
				iEnd += 1
			val possibleLights = 
				for {
					j <- iStart until iEnd
					if targetMatch(xl, fSorted(j))
				} yield {
					val light = fSorted(j)
					var k = j+1
					val heavys = new ArrayBuffer[DinoFeature]
					val hm = light.mass + params.xlinkerDeltaMass
					while (k < fSorted.length &&
							fSorted(k).mass < hm + params.featMzTolerance) {
						val heavy = fSorted(k)
						if (
								light.z == heavy.z && 
								math.abs(light.rtApex - heavy.rtApex - params.xlRtDiff) < params.xlRtDiff95 &&
								math.abs(fSorted(k).mass - hm) < params.featMzTolerance
						)
							heavys += heavy
						k += 1
					}
					(j, light, 
						if (heavys.isEmpty) None
						else Some(heavys.minBy(h => math.abs(light.mass - h.mass)))
					)
				}
			
			XLinkPairs(xl, possibleLights.filter(_._3.isDefined).map(t => FeaturePair(t._1, t._2, t._3.get)))
		}
	}
	
	
	
	def writeMatches(
		results:Seq[XLinkPairs],
		f:File
	)(implicit params:TaxlinkParams) = {
		val w = new BufferedWriter(new FileWriter(f))
		
		def writeRow(qoute:Boolean)(a:Any*) = 
			w.write(a.map(_ match {
				case s:String => 
					if (qoute) params.outQuote + s + params.outQuote
					else s
				case x => x.toString
			}).mkString(params.outSep) + "\n")
			
		
		writeRow(false)(
				"lightMonoMz",
				"heavyMonoMz",
				"rtApexLight",
				"rtApexHeavy",
				"massLight",
				"massHeavy",
				"intLight",
				"intHeavy",
				"intScore",
				"rtScore",
				"massScore",
				"scoreSum",
				"z",
				"xlink"
			)
			
		for {
			XLinkPairs(xl, pairs) <- results
			FeaturePair(i, light, heavy) <- pairs
		} {
			val _intScore = singleIntScore(light.intensitySum, heavy.intensitySum)
			val _rtScore = rtScore(light.rtApex - heavy.rtApex, params.xlRtDiff)
			val _massScore = massScore(light.mass + params.xlinkerDeltaMass, heavy.mass)
			
			writeRow(false)(
					light.mz,
					heavy.mz,
					light.rtApex,
					heavy.rtApex,
					light.mass,
					heavy.mass,
					light.intensityApex,
					heavy.intensityApex,
					_intScore,
					_rtScore,
					_massScore,
					_intScore + _rtScore + _massScore,
					light.z,
					xl.str)
		}
		
		w.close()
	}
	
	
	def singleIntScore(light:Double, heavy:Double)(implicit params:TaxlinkParams) = 
		1 / (1 + math.abs(math.log(light) - math.log(heavy))/params.xlPairInt95)
	
	
	
	def rtScore(target:Double, measured:Double)(implicit params:TaxlinkParams) =
		1.0 / (1 + math.abs(target - measured) / params.xlRtDiff95)
	
	
	
	def massScore(target:Double, measured:Double)(implicit params:TaxlinkParams) =
		1.0 / (1 + math.abs(target - measured) / params.mzAccuracy95)
	
}