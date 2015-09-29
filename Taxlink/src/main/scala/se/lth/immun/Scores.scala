package se.lth.immun

class Scores(
		val mzScore:Double,
		val lightMzScore:Double,
		val heavyMzScore:Double,
		val lightRtScore:Double,
		val heavyRtScore:Double,
		val singleIntScore:Double,
		val doubleIntScore:Double
) {
	override def toString = 
		"mzScore:%.3f   lightMzScore:%.3f   heavyMzScore:%.3f   rtLight:%.3f   rtHeavy:%.3f   intSingle:%.3f   intDouble:%.3f".format(
				mzScore, lightMzScore, heavyMzScore, lightRtScore, heavyRtScore, singleIntScore, doubleIntScore)
}