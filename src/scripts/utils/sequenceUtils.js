import { AMINO_ACIDS, CODONS } from '../configs';


export function initializePixiTexts(textOptions, HGC){

  const codonTexts = {};
  const sequences = Object.keys(AMINO_ACIDS).concat(Object.keys(CODONS));

  sequences.forEach((sequence) => {
    const codonText = {};
    codonText["sequence"] = sequence;

    const pixiText = new HGC.libraries.PIXI.Text(
      sequence,
      textOptions
    );
    pixiText.updateText();

    // We get sharper edges if we scale down a large text
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

export function getAminoAcidsForTile(seq, tileOffset, exonStarts, exonEnds, minX, frontExcessBases){
  const codons = [];
  //let seqFiltered = seq.filter(nuc => nuc !== ".");
  // We cut off the tile offset and last bases, so that we get a length that is a multiple of 3
  //seqFiltered = seqFiltered.slice(tileOffset);
  //seqFiltered = seqFiltered.slice(0, seqFiltered.length - (seqFiltered.length % 3));

  console.log(seq.join(''));
  //console.log(seqFiltered.join(''), seqFiltered.length);
  console.log("tileOffset", tileOffset);

  // Get the first non .
  let firstNonPoint = 0;
  for(let i = 0; i < seq.length; i++){
    if(seq[i] !== "."){
      firstNonPoint = i;
      break;
    }
  }

  // load the first codon
  if(tileOffset > 0){

  }

  let codonSeq = [];
  for(let i = firstNonPoint + tileOffset; i < seq.length; i++){
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
      const codonStr = codonSeq[0].letter + codonSeq[1].letter + codonSeq[2].letter;
      const aa = CODONS[codonStr];
      if(!aa){
        console.warn("Codon " + codonStr + " does not exist. Position: " + codonSeq[0].pos);
        continue;
      }
      // if they are consecutive
      if( (codonSeq[2].pos - codonSeq[1].pos === 1) && (codonSeq[1].pos - codonSeq[0].pos === 1)){
        const codon = {
          posStart: minX + codonSeq[0].pos,
          posEnd: minX + codonSeq[2].pos,
          aminoAcid: aa
        };
        codons.push(codon);
      }
      // Split after first nucleotide
      else if((codonSeq[1].pos - codonSeq[0].pos > 1)){
        const codon1 = {
          posStart: minX + codonSeq[0].pos,
          posEnd: minX + codonSeq[0].pos,
          aminoAcid: aa
        };
        codons.push(codon1);
        const codon2 = {
          posStart: minX + codonSeq[1].pos,
          posEnd: minX + codonSeq[2].pos,
          aminoAcid: aa
        };
        codons.push(codon2);
      }
      // Split after first nucleotide
      else if((codonSeq[2].pos - codonSeq[1].pos > 1)){
        const codon1 = {
          posStart: minX + codonSeq[0].pos,
          posEnd: minX + codonSeq[1].pos,
          aminoAcid: aa
        };
        codons.push(codon1);
        const codon2 = {
          posStart: minX + codonSeq[2].pos,
          posEnd: minX + codonSeq[2].pos,
          aminoAcid: aa
        };
        codons.push(codon2);
      }

      codonSeq = [];
    }
    
  }

  // for(let i = tileOffset; i < seqFiltered.length; i=i+3){
  //   codonStr = seqFiltered[i]+seqFiltered[i+1]+seqFiltered[i+2]
  //   const codon = {
  //     posStart: minX + i
  //   };

  // }

  return codons;

}



export default getContainingExon;
