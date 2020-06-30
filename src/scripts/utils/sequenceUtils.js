import { AMINO_ACIDS, CODONS } from '../configs';


export function initializePixiTexts(textOptions, HGC){

  const codonTexts = {};
  const sequences = Object.keys(AMINO_ACIDS).concat(Object.keys(CODONS));

  Object.keys(CODONS).forEach((sequence) => {
    const codonText = {};
    codonText["sequence"] = sequence;
   
    console.log(sequence, CODONS[sequence].nameAbbrev);
    const pixiText = new HGC.libraries.PIXI.Text(
      CODONS[sequence].nameAbbrev,
      textOptions
    );
    pixiText.updateText();
    
    // We get sharper edges if we scale down a large text
    // This holds the 3 letter AA
    codonText["width"] = pixiText.getBounds().width / 2;
    codonText["height"] = pixiText.getBounds().height / 2;
    codonText["texture"] = pixiText.texture;

    const pixiTextAbbrev = new HGC.libraries.PIXI.Text(
      CODONS[sequence].key,
      textOptions
    );
    pixiTextAbbrev.updateText();

    // We get sharper edges if we scale down a large text
    // This holds the 3 letter AA
    codonText["widthAbbrev"] = pixiTextAbbrev.getBounds().width / 2;
    codonText["heightAbbrev"] = pixiTextAbbrev.getBounds().height / 2;
    codonText["textureAbbrev"] = pixiTextAbbrev.texture;
    

    codonTexts[sequence] = codonText;
  });
  return codonTexts;
}

export function initializePixiTextsOLD(textOptions, HGC){

  const codonTexts = {};
  const sequences = Object.keys(AMINO_ACIDS).concat(Object.keys(CODONS));

  sequences.forEach((sequence) => {
    const codonText = {};
    codonText["sequence"] = sequence;
    let aaText = "";
    if(CODONS[sequence]){
      aaText = CODONS[sequence].nameAbbrev;
    }else if(AMINO_ACIDS[sequence]){
      aaText = sequence;
    }
    console.log(sequence, aaText);
    const pixiText = new HGC.libraries.PIXI.Text(
      aaText,
      textOptions
    );
    pixiText.updateText();
    
    // We get sharper edges if we scale down a large text
    // This holds the 3 letter AA
    codonText["width"] = pixiText.getBounds().width / 2;
    codonText["height"] = pixiText.getBounds().height / 2;
    codonText["texture"] = pixiText.texture;
    

    codonTexts[sequence] = codonText;
  });
  return codonTexts;
}

// Get the exon number, where pos falls in (in chr coordinates)
export function getContainingExon(starts, ends, pos){

  for(let numExon = 0; numExon < starts.length; numExon++){
    if(pos >= starts[numExon] && pos < ends[numExon]){
      const result = {
        exon: numExon,
        start: starts[numExon],
        end: ends[numExon]
      }
      return result;
    }

  }
  return null;
}

export function getTileSequenceOffset(exonStart, exonOffset, tilePos){
  //console.log(exonStart, exonOffset, tilePos);
  const codonStart = exonStart + exonOffset;
  const offset = (3 - ((tilePos - codonStart) % 3)) % 3;
  return offset;
}

export function exonIntersect(starts, ends, startCodonPos, stopCodonPos, evalStart, evalArr){
  const intersection = [...evalArr];
  //console.log("exonIntersect", starts, ends, evalStart, evalArr, intersection.length, evalStart+intersection.length);

  //return [];

  for(let i = 0; i < intersection.length; i++){

    let found = false;
    for(let numExon = 0; numExon < starts.length; numExon++){
      const absPos = i+evalStart;
      if(absPos >= startCodonPos && absPos >= starts[numExon] && absPos < ends[numExon] && absPos < stopCodonPos){
        found = true;
      }
      //console.log(i,starts[numExon],ends[numExon],found);
    }

    if(!found){
      intersection[i] = ".";
    }

  }
  return intersection;
}

export function getNextExon(starts, ends, pos){

  for(let numExon = 0; numExon < starts.length; numExon++){
    if(pos < starts[numExon]){
      const result = {
        exon: numExon,
        start: starts[numExon],
        end: ends[numExon]
      }
      return result;
    }

  }
  return null;
}

export function getAminoAcidsForTile(HGC, seq, tileOffset, chromName, exonStarts, exonEnds, minX, frontExcessBases, pixiTexts, sequenceLoader){
  //console.log("PIXITEXTS", pixiTexts);
  const codons = [];
  let seqFiltered = seq.filter(nuc => nuc !== ".");
  // We cut off the tile offset and last bases, so that we get a length that is a multiple of 3
  //seqFiltered = seqFiltered.slice(tileOffset);
  //seqFiltered = seqFiltered.slice(0, seqFiltered.length - (seqFiltered.length % 3));

  //console.log(seq.join(''));
  //console.log(seqFiltered.join(''), seqFiltered.length);
  //console.log("tileOffset", tileOffset);

  // Get the first non .
  let firstNonPoint = 0;
  for(let i = 0; i < seq.length; i++){
    if(seq[i] !== "."){
      firstNonPoint = i;
      break;
    }
  }

  // Get the last non .
  let lastNonPoint = 0;
  for(let i = seq.length-1; i > 0; i--){
    if(seq[i] !== "."){
      lastNonPoint = i;
      break;
    }
  }

  //console.log("firstNonPoint", firstNonPoint, "lastNonPoint", lastNonPoint)

  let reqNucLeftStart = null;
  let reqNucLeftEnd = null;
  // load the first codon
  if(tileOffset > 0){
    // The required nucleotides are in a previous exon
    if(firstNonPoint > 0){
      const currentExonNum =  getContainingExon(exonStarts, exonEnds, minX + firstNonPoint).exon;
      const previousExon = currentExonNum - 1;
      const previousExonEnd = exonEnds[previousExon];
      reqNucLeftStart = previousExonEnd-(3-tileOffset);
      reqNucLeftEnd = previousExonEnd;
      //requiredNucleotidesStartPromise = sequenceLoader.getSequence(chromName, previousExonEnd-(3-tileOffset), previousExonEnd);
    }
    // The required nucleotides are in the same exon
    else{
      reqNucLeftStart = minX-(3-tileOffset);
      reqNucLeftEnd = minX;
      //requiredNucleotidesStartPromise = sequenceLoader.getSequence(chromName, minX-(3-tileOffset), minX);
      //console.log("ADDITIONSEQ ", requiredNucleotides, minX, firstNonPoint, minX-(3-tileOffset), minX);
    }
    
  }

  let reqNucRightStart = null;
  let reqNucRightEnd = null;
  // load the last codon
  if( (seqFiltered.length + 3 - tileOffset) % 3 !== 0){

    // how many we still need to get
    const numToGet = 3 - ((seqFiltered.length - tileOffset) % 3);

    // We have to get the required nucleotides from the same exon
    if(lastNonPoint === seq.length-1){
      reqNucRightStart = minX+seq.length;
      reqNucRightEnd = minX+seq.length+numToGet;
      //requiredNucleotidesEndPromise = sequenceLoader.getSequence(chromName, minX+seq.length, minX+seq.length+numToGet);
      //console.log("ADDITIONSEQEND ", requiredNucleotides, minX+seq.length, minX+seq.length+numToGet);
    }
    // We have to get the required nucleotides from the next exon
    else{
      const currentExonNum =  getContainingExon(exonStarts, exonEnds, minX + lastNonPoint).exon;
      const nextExon = currentExonNum + 1;
      const nextExonStart = exonStarts[nextExon];
      reqNucRightStart = nextExonStart;
      reqNucRightEnd = nextExonStart+numToGet;
      //console.log("ADDITIONSEQEND2 ", currentExonNum, nextExon, nextExonStart, reqNucRightStart, reqNucRightEnd);
      //requiredNucleotidesEndPromise = sequenceLoader.getSequence(chromName, nextExonStart, nextExonStart+numToGet);
    }
    
  }

  const excessNucleotides = 
  sequenceLoader.getExcessNucleotides(chromName, reqNucLeftStart, reqNucLeftEnd, reqNucRightStart, reqNucRightEnd)
    .then((nucleotides) => {
      let seqMod = seq;
      const initialCodonSeq = [];
      const finalCodonSeq = [];

      for(let i=0; i<nucleotides.length; i++){
        if(nucleotides[i].leftOrRight === "left"){

          const leftNucleotides = nucleotides[i].value.split('');
          //console.log("leftNucleotides", leftNucleotides);
          for(let j=0; j<leftNucleotides.length; j++){
            initialCodonSeq.push({
              pos: reqNucLeftStart - minX + j,
              letter: leftNucleotides[j]
            });
          }
          
        }
        else if(nucleotides[i].leftOrRight === "right"){
          const rightNucleotides = nucleotides[i].value.split('');
          //console.log("rightNucleotides", rightNucleotides);
          for(let j=0; j<rightNucleotides.length; j++){
            finalCodonSeq.push({
              pos: reqNucRightStart - minX + j,
              letter: rightNucleotides[j]
            });
          }

        }
        // Just add the additional Nuleotides to the right of the seqence. getFormattedCodons
        //will take care of it
        // else if(nucleotides[i].leftOrRight === "right"){
        //   seqMod = seqMod.concat(nucleotides[i].value.split(''));
        // }
      }
      //console.log(seq);
      //console.log(reqNucLeftStart, minX);
      //console.log("initialCodonSeq", initialCodonSeq);
      //console.log("finalCodonSeq", finalCodonSeq);
      //console.log("RETURNeED PROMISE ",nucleotides);

      //const codons = getFormattedCodons([], finalCodonSeq, firstNonPoint + tileOffset, seq, pixiTexts, HGC, minX  );
      const codons = getFormattedCodons(initialCodonSeq, finalCodonSeq, firstNonPoint, seq, pixiTexts, HGC, minX  );
      //console.log(codons);
      return codons;
    });

  return excessNucleotides;

}

function getFormattedCodons(initialCodonSeq, finalCodonSeq, seqStart, seq, pixiTexts, HGC, minX ){

  let codons = [];
  let codonSeq = initialCodonSeq;
  for(let i = seqStart; i < seq.length; i++){
    // We are fillig up a codon sequence array. When the length reaches 3, register it and reset.
    if(seq[i] === "."){
      continue;
    }
    //console.log(seq[i]);
    codonSeq.push({
      pos: i,
      letter: seq[i]
    });

    if(codonSeq.length === 3){
      codons = codons.concat(formatCodonSeqence(codonSeq, minX, pixiTexts, HGC));
      codonSeq = [];
    }
    
  }
  //console.log("codonSeq",codonSeq, "finalCodonSeq", finalCodonSeq)
  if(finalCodonSeq.length > 0 &&  codonSeq.length > 0){
    codonSeq = codonSeq.concat(finalCodonSeq);
    if(codonSeq.length === 3){
      codons = codons.concat(formatCodonSeqence(codonSeq, minX, pixiTexts, HGC));
    }
  }

  return codons;
}


function formatCodonSeqence(codonSeq, minX, pixiTexts, HGC){

  const codons = []
  const codonStr = codonSeq[0].letter.toUpperCase() + codonSeq[1].letter.toUpperCase() + codonSeq[2].letter.toUpperCase();
  const aa = CODONS[codonStr];
  if(!aa){
    console.warn("Codon " + codonStr + " does not exist. Position: " + codonSeq[0].pos);
    return codons;
  }

  

  // if they are consecutive
  if( (codonSeq[2].pos - codonSeq[1].pos === 1) && (codonSeq[1].pos - codonSeq[0].pos === 1)){
    const pixiSprite = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].texture);
    pixiSprite.width = pixiTexts[codonStr].width;
    pixiSprite.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].textureAbbrev);
    pixiSpriteAbbrev.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev.height = pixiTexts[codonStr].heightAbbrev;

    const codon = {
      posStart: minX + codonSeq[0].pos,
      posEnd: minX + codonSeq[2].pos,
      aminoAcid: aa,
      width: pixiSprite.width,
      height: pixiSprite.height,
      sprite: pixiSprite,
      widthAbbrev: pixiSpriteAbbrev.width,
      heightAbbrev: pixiSpriteAbbrev.height,
      spriteAbbrev: pixiSpriteAbbrev,
    };
    codons.push(codon);
  }
  // Split after first nucleotide
  else if((codonSeq[1].pos - codonSeq[0].pos > 1)){
    const pixiSprite1 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].texture);
    pixiSprite1.width = pixiTexts[codonStr].width;
    pixiSprite1.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev1 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].textureAbbrev);
    pixiSpriteAbbrev1.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev1.height = pixiTexts[codonStr].heightAbbrev;

    const codon1 = {
      posStart: minX + codonSeq[0].pos,
      posEnd: minX + codonSeq[0].pos,
      aminoAcid: aa,
      width: pixiSpriteAbbrev1.width,
      height: pixiSpriteAbbrev1.height,
      sprite: pixiSpriteAbbrev1,
      widthAbbrev: pixiSpriteAbbrev1.width,
      heightAbbrev: pixiSpriteAbbrev1.height,
      spriteAbbrev: pixiSpriteAbbrev1,
    };
    codons.push(codon1);

    const pixiSprite2 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].texture);
    pixiSprite2.width = pixiTexts[codonStr].width;
    pixiSprite2.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev2 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].textureAbbrev);
    pixiSpriteAbbrev2.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev2.height = pixiTexts[codonStr].heightAbbrev;

    const codon2 = {
      posStart: minX + codonSeq[1].pos,
      posEnd: minX + codonSeq[2].pos,
      aminoAcid: aa,
      width: pixiSprite2.width,
      height: pixiSprite2.height,
      sprite: pixiSprite2,
      widthAbbrev: pixiSpriteAbbrev2.width,
      heightAbbrev: pixiSpriteAbbrev2.height,
      spriteAbbrev: pixiSpriteAbbrev2,
    };
    codons.push(codon2);
  }
  // Split after first nucleotide
  else if((codonSeq[2].pos - codonSeq[1].pos > 1)){
    const pixiSprite1 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].texture);
    pixiSprite1.width = pixiTexts[codonStr].width;
    pixiSprite1.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev1 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].textureAbbrev);
    pixiSpriteAbbrev1.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev1.height = pixiTexts[codonStr].heightAbbrev;

    const codon1 = {
      posStart: minX + codonSeq[0].pos,
      posEnd: minX + codonSeq[1].pos,
      aminoAcid: aa,
      width: pixiSprite1.width, 
      height: pixiSprite1.height,
      sprite: pixiSprite1,
      widthAbbrev: pixiSpriteAbbrev1.width,
      heightAbbrev: pixiSpriteAbbrev1.height,
      spriteAbbrev: pixiSpriteAbbrev1,
    };
    codons.push(codon1);

    const pixiSprite2 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].texture);
    pixiSprite2.width = pixiTexts[codonStr].width;
    pixiSprite2.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev2 = new HGC.libraries.PIXI.Sprite(pixiTexts[codonStr].textureAbbrev);
    pixiSpriteAbbrev2.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev2.height = pixiTexts[codonStr].heightAbbrev;

    const codon2 = {
      posStart: minX + codonSeq[2].pos,
      posEnd: minX + codonSeq[2].pos,
      aminoAcid: aa,
      width: pixiSpriteAbbrev2.width, // Display the short version by default
      height: pixiSpriteAbbrev2.height,
      sprite: pixiSpriteAbbrev2,
      widthAbbrev: pixiSpriteAbbrev2.width,
      heightAbbrev: pixiSpriteAbbrev2.height,
      spriteAbbrev: pixiSpriteAbbrev2,
    };
    codons.push(codon2);
  }

  return codons;


}


export default getContainingExon;
