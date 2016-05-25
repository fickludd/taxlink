package se.lth.immun

import java.util.Properties
import java.io.File
import se.jt.CLIApp
import se.lth.immun.xlink.XLink
import se.lth.immun.unimod.UniMod
import se.lth.immun.chem._

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

object XLAssayGen extends CLIApp {
	
	
	val C13C12_DIFF = 1.0033548378
	//val BY_FRAGS = Array(EPeptideFragment.b, EPeptideFragment.y)
	
	abstract class XLPeptidePos(val top:Boolean, val nTerm:Boolean)
	case object TopNTerm extends XLPeptidePos(true, true)
	case object TopCTerm extends XLPeptidePos(true, false)
	case object BottomNTerm extends XLPeptidePos(false, true)
	case object BottomCTerm extends XLPeptidePos(false, false)
	
	case class PrecursorIsotopeTarget(mz:Double, order:Int, occurence:Double, label:String)
	case class FragmentTransition(
			precMz:Double, fragMz:Double, 
			fragmentType:Char, ordinal:Int, pos:Int,
			peptide:XLPeptidePos)
	
	case class XLAssay(
			xlinkStr:String,
			xlink:XLink,
			z:Int,
			mz:Double,
			rtLight:Double,
			txIntScore:Double,
			txRtScore:Double,
			txMassScore:Double,
			kjQValue:Double,
			targets:Seq[PrecursorIsotopeTarget],
			transitions:Seq[FragmentTransition]
			)

	def main(args:Array[String]):Unit = {
		var properties = new Properties
	    properties.load(this.getClass.getResourceAsStream("/pom.properties"))
		val appName 		= properties.getProperty("pom.artifactId")
		val appVersion 	= properties.getProperty("pom.version")
    
		implicit val params = new XLAssayGenParams
		
		failOnError(parseArgs(appName, appVersion, args, params, List("xlMaster"), None))
		
		params.parseFragTypes
		
		println("   "+appName+" "+appVersion)
		println("   xl master file: " + params.xlMaster.value)
		println("   fragment types: " + params.FRAG_TYPES.mkString(","))
		
		val masterRecords = XLMasterFormat.read(params.xlMaster)
		
		val filteredRecords = masterRecords.filter(r => 
			r.kojak.qValue < params.qValue && (r.isComplete || params.useIncompleteRows))
		
		val xlAssayTries = filteredRecords.map(genAssays(params))
		val xlAssays = xlAssayTries.collect({
				case Left(assay) 	=> Some(assay)
				case Right(err) 	=> println(err); None
			}).flatten
		
		val (dir, outName) = params.outBase
		if (params.outHumanImg)
			HumanImg.report(new File(dir, outName), xlAssays.take(10))
			
		XLTraml.write(new File(dir, outName + ".traml"), xlAssays)
	}
	
	
	
	def genAssays(params:XLAssayGenParams)(r:XLMasterFormat.Row):Either[XLAssay, String] = {
		XLink.fromString(r.common.xlink, UniMod.parseUniModSequence, params.kojakConf) match {
			case Left(xlink) =>
				val targets = genTargets(xlink, r.common.z, params)
				Left(XLAssay(
						r.common.xlink, 
						xlink, 
						r.common.z, 
						r.common.mz,
						r.taxlink.rtApexLight,
						r.taxlink.intScore,
						r.taxlink.rtScore,
						r.taxlink.massScore,
						r.kojak.qValue,
						targets,
						genTransitions(xlink, targets.head, r.common.z, params)
					))
				
			case Right(err) => Right(err)
				
		}
	}
	
	
	def genTargets(xlink:XLink, z:Int, params:XLAssayGenParams) = {
		val lightXL = xlink.getComposition
		val targets = new ArrayBuffer[PrecursorIsotopeTarget]
		
		def addIsotopesUntil(atoms:ElementComposition, label:String) = {
			val m0 = atoms.monoisotopicMass
			val inOrderOfOccurence = new Queue[(Double, Int)]
			inOrderOfOccurence ++= lightXL.getIsotopeDistribution.intensities.zipWithIndex.sortBy(_._1).reverse
			var occurenceMass = 0.0
			while (occurenceMass < params.isotopeOccurenceMass) {
				val (occ, ind) = inOrderOfOccurence.dequeue
				val order = ind
				val m = m0 + C13C12_DIFF * order
				occurenceMass += occ
				targets += PrecursorIsotopeTarget((m/z + Constants.PROTON_WEIGHT), order, occ, label)
			}
		}
		
		addIsotopesUntil(lightXL, "d0")
		if (params.xlDeuteriums > 0) {
			import Element._
			val n:Int = params.xlDeuteriums
			val heavyLinkerDiff = new ElementComposition(Array(H, H2), Array(-n, n))
			val heavyXL = lightXL.join(heavyLinkerDiff)
			addIsotopesUntil(heavyXL, "d"+params.xlDeuteriums.value)
		}
		
		targets
	}
	
	
	def genTransitions(
			xlink:XLink, 
			precTarget:PrecursorIsotopeTarget,
			z:Int, 
			params:XLAssayGenParams
	) = {
		val precMass = xlink.monoisotopicMass + C13C12_DIFF * precTarget.order
		val precMz = m2mz(precMass, z)
		
		val transitions = new ArrayBuffer[FragmentTransition]
		
		val positions = xlink.positions.flatMap(link => 
				List(link.pepPos1, link.pepPos2)
			).collect { 
				case kp:XLink.KnownPosition => kp 
			}
		
		def toFragTrans(xlPepPos:XLPeptidePos, peptideLength:Int)(pf:PeptideFragment) =
			FragmentTransition(
				precMz, 
				m2mz(pf.mass, 1), 
				pf.fragmentType.name.head,
				pf.ordinal,
				if (pf.fragmentType == EPeptideFragment.b) pf.ordinal
				else peptideLength - pf.ordinal,
				xlPepPos
			)
			
		val (topN, topC) = genPepFragments(xlink.pep1, positions, params.FRAG_TYPES)
		transitions ++= topN.map(toFragTrans(TopNTerm, xlink.pep1.aminoAcids.length))
		transitions ++= topC.map(toFragTrans(TopCTerm, xlink.pep1.aminoAcids.length))
		
		for (pep2 <- xlink.pep2) {
			val (botN, botC) = genPepFragments(pep2, positions, params.FRAG_TYPES)
			transitions ++= botN.map(toFragTrans(BottomNTerm, pep2.aminoAcids.length))
			transitions ++= botC.map(toFragTrans(BottomNTerm, pep2.aminoAcids.length))
		}
		
		transitions
	}
	
	
	def genPepFragments(
			pep:Peptide,
			positions:Seq[XLink.KnownPosition],
			fragTypes:Array[EPeptideFragment]
	) = {
		import EPeptideFragment._
		
		val pos = positions.filter(_.pep == pep).map(_.pos)
		val frags = pep.getFragments(fragTypes)
		val nTermFrags = frags.filter(pf => pf.fragmentType.isNTerm && pf.ordinal < pos.min)
		val cTermFrags = frags.filter(pf => !pf.fragmentType.isNTerm && pos.max <= pep.aminoAcids.length - pf.ordinal)
		(nTermFrags, cTermFrags)
	}
	
	
	def m2mz(m:Double, z:Int) = m / z + Constants.PROTON_WEIGHT
}