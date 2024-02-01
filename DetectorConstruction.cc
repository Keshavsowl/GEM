//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Tubs.hh"
// #include "G4MaterialMixture.hh"



namespace B1
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::DetectorConstruction()
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::~DetectorConstruction()
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::Construct()  
  {
    // Get nist material manager
    G4cout << "Defining materials" << G4endl;
    G4NistManager *nist = G4NistManager::Instance();
    
     G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");
     G4Material *box_Cu = nist->FindOrBuildMaterial("G4_Cu");
     G4Material *box_Ar = nist->FindOrBuildMaterial("G4_Ar");
     G4Material *box_Kapton = nist->FindOrBuildMaterial("G4_Kapton");
    
     
     
  /*   G4Element* elAr = nistManager->FindOrBuildElement("Ar");
     G4Element* elC = nistManager->FindOrBuildElement("C");
     G4Element* elO = nistManager->FindOrBuildElement("O");

     // Define the material for the mixture of CO2 and Ar
     G4double density = 1.782e-3 * g/cm3;  // Density of the mixture
     G4Material* ArCO2Mixture = new G4Material("ArCO2Mixture", density, 2);
     ArCO2Mixture->AddMaterial(new G4Material("CO2", density, 2), 0.7);          // 70% CO2
     ArCO2Mixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ar"), 0.3);  // 30% Argon
     ArCO2Mixture->GetElement(0)->AddIsotope(elC, 1);  // One atom of Carbon
     ArCO2Mixture->GetElement(0)->AddIsotope(elO, 2);  // Two atoms of Oxygen  */
     
     
  
  
       
  
  
  
  
     

   G4Material *box_kapton = nist->FindOrBuildMaterial("G4_Kapton");
    
    // Option to switch on/off checking of volumes overlaps
        G4cout << "Defining world volume" << G4endl;
    G4bool checkOverlaps = true;

    // construct a world volume box of dimensions 100 cm x 100 cm x 100 cm
    G4Box  *solidWorld = new G4Box("World", 100* cm, 100*cm, 100*cm);                    // Note convention: Half length is used
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World"); // its name
        G4cout << "Defining PVPworld" << G4endl;
    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (0,0,0)
                          logicWorld,      // its logical volume
                          "World",         // its name
                          0,               // its mother  volume
                          false,           // no boolean operation
                          0,               // copy number
                          checkOverlaps);  // overlaps checking
                          
                          
                          
     
      
 
    G4Box* outerBox = new G4Box("OuterBox", 50 * cm, 50 * cm, 50 * cm);
    G4Box* innerBox = new G4Box("InnerBox", 48 * cm, 48 * cm, 48 * cm);
    G4SubtractionSolid* hollowBox = new G4SubtractionSolid("HollowBox", outerBox, innerBox);
    G4LogicalVolume* logicvolume0 = new G4LogicalVolume(hollowBox, box_Ar , "GasVolume");
    G4VPhysicalVolume* physGas = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, 0),
                                               logicvolume0,
                                               "GasVolume",
                                               logicWorld,
                                               false,
                                               0,
                                               checkOverlaps);

         
      // Translation and Rotation Vector  
      
    G4ThreeVector zzTrans(0, 0, 0);
    
    G4RotationMatrix* yRot =  new G4RotationMatrix;
    yRot -> rotateY(0*deg);                       
    
    
    
   
    
    // Defining Layer of Kapton
                             
    G4Box *box = new G4Box ("Box", 20*cm, 20*cm, 25* um);                   
    G4Tubs* cylinderSolid = new G4Tubs("Cylinder", 0*cm , 2*cm, 20*cm, 0*degree, 360* degree);
    G4SubtractionSolid* subtraction1 = new G4SubtractionSolid("Box-Cylinder", box , cylinderSolid, yRot, zzTrans);
    for ( int i= -4; i<=4; i++ ){
    
    for (int j = -4; j<=4; j++) {
    
    G4ThreeVector zzTrans(4*j*cm,4*i*cm, 0);
    
    subtraction1 = new G4SubtractionSolid("Box-Cylinder", subtraction1 , cylinderSolid, yRot, zzTrans);
    } }
      
      G4LogicalVolume* logicBox1 =   new G4LogicalVolume(subtraction1, box_Kapton,"Box-box"); 
      G4VPhysicalVolume *physBox1 =
      new G4PVPlacement(0,                      // rotation matrix
                          G4ThreeVector(0, 0, 0 ), // at // translation
                          logicBox1,               // its logical volume
                          "Box-Cube2",         // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operation
                          0,                      // copy number
                          checkOverlaps);         // overlaps checking 
                          
                          
   
   
   
   
   
   // Defining layer of Copper
    
   G4Box *box2 = new G4Box ("Box", 20*cm, 20*cm, 2.5* um);                   
    G4SubtractionSolid* subtraction2 = new G4SubtractionSolid("Box-Cylinder", box2 , cylinderSolid, yRot, zzTrans);
    for ( int i=1; i<=1 ; i++ ){
    
    for (int j =1; j<=1; j++) {
    
    G4ThreeVector zzTrans(4*j*cm , 4*i*cm, 0);
    
    subtraction2 = new G4SubtractionSolid("Box-Cylinder", subtraction2 , cylinderSolid, yRot, zzTrans);
    } }
       
  G4LogicalVolume* logicBox2 =   new G4LogicalVolume(subtraction1,box_Cu,"Box-box"); 
   G4VPhysicalVolume *physBox2 =
                          new G4PVPlacement(0,                      // rotation matrix
                          G4ThreeVector(0*cm, 0*cm, 27.5*cm ), // at // translation
                          logicBox2,               // its logical volume
                          "Box-Cube2",         // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operation
                          0,                      // copy number
                          checkOverlaps);         // overlaps checking 
                          
                       
    G4VPhysicalVolume *physBox3 =
        new G4PVPlacement(0,                      // rotation matrix
                          G4ThreeVector(0*cm, 0*cm, -27.5*cm), // at // translation
                          logicBox2,               // its logical volume
                          "Box-Cube3",         // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operation
                          0,                      // copy number
                          checkOverlaps);         // overlaps checking         */                 
                          
               
    
  //  fScoringVolume = logicBox;

    //
    // always return the physical World
    //
    return physWorld;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
