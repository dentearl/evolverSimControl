import unittest
import libSimControl as LSC
import libSimTree as LST

class VerifyProgramsBadInput( unittest.TestCase ):
    def test_InputNotList( self ):
        """input to verifyPrograms must a list"""
        self.assertRaises( LSC.BadInputError, LSC.verifyPrograms, 'banana_432215234231113')
        self.assertRaises( LSC.BadInputError, LSC.verifyPrograms, 50)
    def test_InputListContainsNotStrings( self ):
        """All members of the input list must be strings"""
        self.assertRaises( LSC.BadInputError, LSC.verifyPrograms, [ 20, 5.0])
    def test_ProgramDoesNotExist( self ):
        """verifyPrograms should fail if a program does not exist"""
        self.assertRaises( LSC.ProgramDoesNotExistError, LSC.verifyPrograms,
                           ['_aoeu_AOEU_SNTH_snth_','I dont think this would exist __1_'] )

class CheckCommandPacker( unittest.TestCase ):
    knownValues = (('echo "$SHELL"', '&myCMD;echo &myQuot;$SHELL&myQuot;&myCMD;'),
                   ('#echo "horses" // \'\'', '&myCMD;#echo &myQuot;horses&myQuot; // &myApos;&myApos;&myCMD;'))
    def test_BadInput( self ):
        """commandPacker takes a single string as input"""
        self.assertRaises( LSC.BadInputError, LSC.commandPacker, ['hi there'])
        self.assertRaises( LSC.BadInputError, LSC.commandPacker, 20.0 )
        self.assertRaises( LSC.BadInputError, LSC.commandPacker, {'steve':'miranda'} )
    def test_commandPackerKnownValues( self ):
        """commandPacker should give known result with known input"""
        for input, output in self.knownValues:
            self.assertEqual( output, LSC.commandPacker( input ) )

class CheckFixName( unittest.TestCase ):
    knownValues = (('((simGorilla:0.008825,(simHuman:0.006700,simChimp:0.006667)sHuman-sChimp:0.002250)sG-sH-sC:0.009680,simOrang:0.018318):0.000000;',
                    '_L__L_simGorilla-0_008825_L_simHuman-0_006700simChimp-0_006667_R_sHuman-sChimp-0_002250_R_sG-sH-sC-0_009680simOrang-0_018318_R'),
                   ('(simHuman:0.002)', '_L_simHuman-0_002_R'),
                   ('',''))
    def test_fixNameKnownValues( self ):
        """fixName should return known values for known input"""
        for input, output in self.knownValues:
            self.assertEqual( output, LST.fixName(input) )

if __name__ == "__main__":
    unittest.main()
