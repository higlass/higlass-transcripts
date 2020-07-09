import { AMINO_ACIDS, CODONS } from "../configs";

export function initializePixiTexts(textOptions, HGC) {
  const codonTexts = {};
  //const sequences = Object.keys(AMINO_ACIDS).concat(Object.keys(CODONS));

  Object.keys(CODONS).forEach((sequence) => {
    const codonText = {};
    codonText["sequence"] = sequence;

    //console.log(sequence, CODONS[sequence].nameAbbrev);
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


// Get the exon number, where pos falls in (in chr coordinates)
export function getContainingExon(starts, ends, pos) {
  for (let numExon = 0; numExon < starts.length; numExon++) {
    if (pos >= starts[numExon] && pos < ends[numExon]) {
      const result = {
        exon: numExon,
        start: starts[numExon],
        end: ends[numExon],
      };
      return result;
    }
  }
  return null;
}

export function getTileSequenceOffset(exonStart, exonOffset, tilePos) {
  //console.log(exonStart, exonOffset, tilePos);
  const codonStart = exonStart + exonOffset;
  const offset = (3 - ((tilePos - codonStart) % 3)) % 3;
  return offset;
}

export function exonIntersect(
  starts,
  ends,
  startCodonPos,
  stopCodonPos,
  evalStart,
  evalArr
) {
  const intersection = [...evalArr];
  
  for (let i = 0; i < intersection.length; i++) {
    let found = false;
    for (let numExon = 0; numExon < starts.length; numExon++) {
      const absPos = i + evalStart;
      if (
        absPos >= startCodonPos &&
        absPos >= starts[numExon] &&
        absPos < ends[numExon] &&
        absPos < stopCodonPos
      ) {
        found = true;
      }
    }

    if (!found) {
      intersection[i] = ".";
    }
  }
  return intersection;
}

export function getNextExon(starts, ends, pos) {
  for (let numExon = 0; numExon < starts.length; numExon++) {
    if (pos < starts[numExon]) {
      const result = {
        exon: numExon,
        start: starts[numExon],
        end: ends[numExon],
      };
      return result;
    }
  }
  return null;
}

export function reverseChrCoord(pos, chromeSize) {
  if (Array.isArray(pos)) {
    const newPos = [];
    for (let i = 0; i < pos.length; i++) {
      newPos.push(chromeSize - pos[i]);
    }
    return newPos.sort();
  } else {
    return chromeSize - pos;
  }
}

export function getMinusStrandSeq(seq) {
  const seqSplit = seq.split("");
  const newSeq = [];
  for (let i = seqSplit.length - 1; i >= 0; i--) {
    if (seqSplit[i].toUpperCase() === "A") {
      newSeq.push("T");
    } else if (seqSplit[i].toUpperCase() === "T") {
      newSeq.push("A");
    } else if (seqSplit[i].toUpperCase() === "C") {
      newSeq.push("G");
    } else if (seqSplit[i].toUpperCase() === "G") {
      newSeq.push("C");
    } else {
      console.warn(
        "Couldn't convert " + seqSplit[i] + " in minus strand conversion"
      );
    }
  }

  return newSeq.join("");
}

export function getAminoAcidsForTile(
  HGC,
  seq,
  tileOffset,
  chromName,
  chromLength,
  strand,
  exonStarts,
  exonEnds,
  minX,
  pixiTexts,
  sequenceLoader
) {

  const codons = [];
  let seqFiltered = seq.filter((nuc) => nuc !== ".");
  // We cut off the tile offset and last bases, so that we get a length that is a multiple of 3
 
  // Get the first non .
  let firstNonPoint = 0;
  for (let i = 0; i < seq.length; i++) {
    if (seq[i] !== ".") {
      firstNonPoint = i;
      break;
    }
  }

  // Get the last non .
  let lastNonPoint = 0;
  for (let i = seq.length - 1; i > 0; i--) {
    if (seq[i] !== ".") {
      lastNonPoint = i;
      break;
    }
  }

  //console.log("firstNonPoint", firstNonPoint, "lastNonPoint", lastNonPoint, minX)

  let reqNucLeftStart = null;
  let reqNucLeftEnd = null;
  // load the first codon
  if (tileOffset > 0) {
    // The required nucleotides are in a previous exon
    if (firstNonPoint > 0) {
      const currentExonNum = getContainingExon(
        exonStarts,
        exonEnds,
        minX + firstNonPoint
      ).exon;

      const previousExon = currentExonNum - 1;
      const previousExonEnd = exonEnds[previousExon];
      reqNucLeftStart = previousExonEnd - (3 - tileOffset);
      reqNucLeftEnd = previousExonEnd;
    }
    // The required nucleotides are in the same exon
    else {
      reqNucLeftStart = minX - (3 - tileOffset);
      reqNucLeftEnd = minX;
    }
  }

  let reqNucRightStart = null;
  let reqNucRightEnd = null;
  // load the last codon
  if ((seqFiltered.length + 3 - tileOffset) % 3 !== 0) {
    // how many we still need to get
    const numToGet = 3 - ((seqFiltered.length - tileOffset) % 3);

    // We have to get the required nucleotides from the same exon
    if (lastNonPoint === seq.length - 1) {
      reqNucRightStart = minX + seq.length;
      reqNucRightEnd = minX + seq.length + numToGet;
    }
    // We have to get the required nucleotides from the next exon
    else {
      const currentExonNum = getContainingExon(
        exonStarts,
        exonEnds,
        minX + lastNonPoint
      ).exon;
      const nextExon = currentExonNum + 1;
      const nextExonStart = exonStarts[nextExon];
      reqNucRightStart = nextExonStart;
      reqNucRightEnd = nextExonStart + numToGet;
    }
  }

  const excessNucleotides = sequenceLoader
    .getExcessNucleotides(
      chromName,
      chromLength,
      strand,
      reqNucLeftStart,
      reqNucLeftEnd,
      reqNucRightStart,
      reqNucRightEnd
    )
    .then((nucleotides) => {

      const initialCodonSeq = [];
      const finalCodonSeq = [];

      for (let i = 0; i < nucleotides.length; i++) {
        if (nucleotides[i].leftOrRight === "left") {
          const leftNucleotides = nucleotides[i].value.split("");

          for (let j = 0; j < leftNucleotides.length; j++) {
            initialCodonSeq.push({
              pos: reqNucLeftStart - minX + j,
              letter: leftNucleotides[j],
            });
          }
        } else if (nucleotides[i].leftOrRight === "right") {
          const rightNucleotides = nucleotides[i].value.split("");

          for (let j = 0; j < rightNucleotides.length; j++) {
            finalCodonSeq.push({
              pos: reqNucRightStart - minX + j,
              letter: rightNucleotides[j],
            });
          }
        }
        
      }
     
      const codons = getFormattedCodons(
        initialCodonSeq,
        finalCodonSeq,
        firstNonPoint,
        seq,
        pixiTexts,
        strand,
        chromLength,
        HGC,
        minX
      );
      return codons;
    });


  return excessNucleotides;
}

function getFormattedCodons(
  initialCodonSeq,
  finalCodonSeq,
  seqStart,
  seq,
  pixiTexts,
  strand,
  chromLength,
  HGC,
  minX
) {
  let codons = [];
  let codonSeq = initialCodonSeq;
  for (let i = seqStart; i < seq.length; i++) {
    // We are fillig up a codon sequence array. When the length reaches 3, register it and reset.
    if (seq[i] === ".") {
      continue;
    }

    codonSeq.push({
      pos: i,
      letter: seq[i],
    });

    if (codonSeq.length === 3) {
      codons = codons.concat(
        formatCodonSeqence(codonSeq, minX, pixiTexts, strand, chromLength, HGC)
      );
      codonSeq = [];
    }
  }

  if (finalCodonSeq.length > 0 && codonSeq.length > 0) {
    codonSeq = codonSeq.concat(finalCodonSeq);
    if (codonSeq.length === 3) {
      codons = codons.concat(
        formatCodonSeqence(codonSeq, minX, pixiTexts, strand, chromLength, HGC)
      );
    }
  }

  return codons;
}

function formatCodonSeqence(
  codonSeq,
  minX,
  pixiTexts,
  strand,
  chromLength,
  HGC
) {
  const codons = [];
  const codonStr =
    codonSeq[0].letter.toUpperCase() +
    codonSeq[1].letter.toUpperCase() +
    codonSeq[2].letter.toUpperCase();
  const aa = CODONS[codonStr];
  if (!aa) {
    console.warn(
      "Codon " + codonStr + " does not exist. Position: " + codonSeq[0].pos
    );
    return codons;
  }

  // if they are consecutive
  if (
    codonSeq[2].pos - codonSeq[1].pos === 1 &&
    codonSeq[1].pos - codonSeq[0].pos === 1
  ) {
    const pixiSprite = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].texture
    );
    pixiSprite.width = pixiTexts[codonStr].width;
    pixiSprite.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].textureAbbrev
    );
    pixiSpriteAbbrev.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev.height = pixiTexts[codonStr].heightAbbrev;
    const posStart =
      strand === "+"
        ? minX + codonSeq[0].pos
        : reverseChrCoord(minX + codonSeq[2].pos + 1, chromLength);
    const posEnd =
      strand === "+"
        ? minX + codonSeq[2].pos
        : reverseChrCoord(minX + codonSeq[0].pos + 1, chromLength);

    const codon = {
      posStart: posStart,
      posEnd: posEnd,
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
  else if (codonSeq[1].pos - codonSeq[0].pos > 1) {
    const pixiSprite1 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].texture
    );
    pixiSprite1.width = pixiTexts[codonStr].width;
    pixiSprite1.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev1 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].textureAbbrev
    );
    pixiSpriteAbbrev1.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev1.height = pixiTexts[codonStr].heightAbbrev;

    const posStart1 =
      strand === "+"
        ? minX + codonSeq[0].pos
        : reverseChrCoord(minX + codonSeq[0].pos + 1, chromLength);
    const posEnd1 =
      strand === "+"
        ? minX + codonSeq[0].pos
        : reverseChrCoord(minX + codonSeq[0].pos + 1, chromLength);

    const codon1 = {
      posStart: posStart1,
      posEnd: posEnd1,
      aminoAcid: aa,
      width: pixiSpriteAbbrev1.width,
      height: pixiSpriteAbbrev1.height,
      sprite: pixiSpriteAbbrev1,
      widthAbbrev: pixiSpriteAbbrev1.width,
      heightAbbrev: pixiSpriteAbbrev1.height,
      spriteAbbrev: pixiSpriteAbbrev1,
    };
    codons.push(codon1);

    const pixiSprite2 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].texture
    );
    pixiSprite2.width = pixiTexts[codonStr].width;
    pixiSprite2.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev2 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].textureAbbrev
    );
    pixiSpriteAbbrev2.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev2.height = pixiTexts[codonStr].heightAbbrev;

    const posStart2 =
      strand === "+"
        ? minX + codonSeq[1].pos
        : reverseChrCoord(minX + codonSeq[2].pos + 1, chromLength);
    const posEnd2 =
      strand === "+"
        ? minX + codonSeq[2].pos
        : reverseChrCoord(minX + codonSeq[1].pos + 1, chromLength);

    const codon2 = {
      posStart: posStart2,
      posEnd: posEnd2,
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
  else if (codonSeq[2].pos - codonSeq[1].pos > 1) {
    const pixiSprite1 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].texture
    );
    pixiSprite1.width = pixiTexts[codonStr].width;
    pixiSprite1.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev1 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].textureAbbrev
    );
    pixiSpriteAbbrev1.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev1.height = pixiTexts[codonStr].heightAbbrev;

    const posStart1 =
      strand === "+"
        ? minX + codonSeq[0].pos
        : reverseChrCoord(minX + codonSeq[1].pos + 1, chromLength);
    const posEnd1 =
      strand === "+"
        ? minX + codonSeq[1].pos
        : reverseChrCoord(minX + codonSeq[0].pos + 1, chromLength);

    const codon1 = {
      posStart: posStart1,
      posEnd: posEnd1,
      aminoAcid: aa,
      width: pixiSprite1.width,
      height: pixiSprite1.height,
      sprite: pixiSprite1,
      widthAbbrev: pixiSpriteAbbrev1.width,
      heightAbbrev: pixiSpriteAbbrev1.height,
      spriteAbbrev: pixiSpriteAbbrev1,
    };
    codons.push(codon1);

    const pixiSprite2 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].texture
    );
    pixiSprite2.width = pixiTexts[codonStr].width;
    pixiSprite2.height = pixiTexts[codonStr].height;

    const pixiSpriteAbbrev2 = new HGC.libraries.PIXI.Sprite(
      pixiTexts[codonStr].textureAbbrev
    );
    pixiSpriteAbbrev2.width = pixiTexts[codonStr].widthAbbrev;
    pixiSpriteAbbrev2.height = pixiTexts[codonStr].heightAbbrev;

    const posStart2 =
      strand === "+"
        ? minX + codonSeq[2].pos
        : reverseChrCoord(minX + codonSeq[2].pos + 1, chromLength);
    const posEnd2 =
      strand === "+"
        ? minX + codonSeq[2].pos
        : reverseChrCoord(minX + codonSeq[2].pos + 1, chromLength);

    const codon2 = {
      posStart: posStart2,
      posEnd: posEnd2,
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
