import { RemoteFile } from "generic-filehandle";
import {
  reverseChrCoord,
  getMinusStrandSeq,
} from "./utils";

class SequenceLoader {
  //constructor(dataConfig) {
  constructor(fastaUrl, faiUrl) {
    //this.dataConfig = dataConfig;
    //this.trackUid = slugid.nice();

    const { IndexedFasta } = require("@gmod/indexedfasta");

    const remoteFA = new RemoteFile(fastaUrl);
    const remoteFAI = new RemoteFile(faiUrl);

    this.sequenceFile = new IndexedFasta({
      fasta: remoteFA,
      fai: remoteFAI,
    });
  }


  // We assume that we are looking for subsequences within a chromosome
  getSubSequence(chromName, exonStarts, exonEnds, startCodonPos, stopCodonPos){

    const recordPromises = [];

    for(let i=0; i < exonStarts.length; i++){
      const curStart = Math.min(Math.max(exonStarts[i], startCodonPos), exonEnds[i]);
      const curEnd = Math.max(Math.min(exonEnds[i], stopCodonPos), exonStarts[i]);
      if(curStart >= curEnd){
        continue;
      }

      recordPromises.push(
        this.sequenceFile
          .getSequence(
            chromName,
            curStart,
            curEnd
          )
          .then((value) => {
            return value;
          })
      );
    }
    

    return Promise.all(recordPromises).then((values) => {
      return values.join('');
    });

  }

  // We assume that we are looking for subsequences within a chromosome
  getExcessNucleotides(chromName, chromLength, strand, start1, end1, start2, end2){

    //console.log(chromName, chromLength, strand, start1, end1, start2, end2)

    let start1used = start1;
    let end1used = end1;
    let start2used = start2;
    let end2used = end2;

    if(strand === "-"){
      start1used = reverseChrCoord(end1, chromLength);
      end1used = reverseChrCoord(start1, chromLength);
      start2used = reverseChrCoord(end2, chromLength);;
      end2used = reverseChrCoord(start2, chromLength);;
    }

    const recordPromises = [];

    if(start1used && end1used){
      recordPromises.push(
        this.sequenceFile
          .getSequence(
            chromName,
            start1used,
            end1used
          )
          .then((value) => {
            return {
              leftOrRight: "left",
              value: strand === "+" ? value : getMinusStrandSeq(value)
            };
          })
      );
    }

    if(start2used && end2used){
      recordPromises.push(
        this.sequenceFile
          .getSequence(
            chromName,
            start2used,
            end2used
          )
          .then((value) => {
            return {
              leftOrRight: "right",
              value: strand === "+" ? value : getMinusStrandSeq(value)
            };
          })
      );
    }

    
    return Promise.all(recordPromises).then((values) => {
      return values;
    });

  }

  // We assume that we are looking for subsequences within a chromosome
  getSequence(chromName, start, end){

    const recordPromises = [];

    recordPromises.push(
      this.sequenceFile
        .getSequence(
          chromName,
          start,
          end
        )
        .then((value) => {
          return value;
        })
    );
    
    return Promise.all(recordPromises).then((values) => {
      return values;
    });

  }

  // get the sequence for a given tile with optionally an additional number of nuleodides in the beginning
  getTile(z, x, tsInfo) {

    const tileWidth = +tsInfo.max_width / 2 ** +z;

    // get the bounds of the tile
    let minX = tsInfo.min_pos[0] + x * tileWidth;
    const maxX = tsInfo.min_pos[0] + (x + 1) * tileWidth;

    const chromSizes = tsInfo.chrom_sizes.split('\t').map(x=>+x);
    const chromNames = tsInfo.chrom_names.split('\t');

    const recordPromises = [];

    this.chromInfo = this.parseChromsizes(chromNames, chromSizes);

    const { chromLengths, cumPositions } = this.chromInfo;

    for (let i = 0; i < cumPositions.length; i++) {
      const chromName = cumPositions[i].chr;
      const chromStart = cumPositions[i].pos;

      const chromEnd = cumPositions[i].pos + chromLengths[chromName];

      if (chromStart <= minX && minX < chromEnd) {
        // start of the visible region is within this chromosome

        if (maxX > chromEnd) {
          // the visible region extends beyond the end of this chromosome
          // fetch from the start until the end of the chromosome
          recordPromises.push(
            this.sequenceFile
              .getSequence(
                chromName,
                minX - chromStart,
                chromEnd - chromStart
              )
              .then((value) => {
                return value;
              })
          );

          // continue onto the next chromosome
          minX = chromEnd;
        } else {
          const endPos = Math.ceil(maxX - chromStart);
          const startPos = Math.floor(minX - chromStart);
          // the end of the region is within this chromosome
          recordPromises.push(
            this.sequenceFile
              .getSequence(chromName, startPos, endPos)
              .then((value) => {
                return value;
              })
          );

          // end the loop because we've retrieved the last chromosome
          break;
        }
      }
    }

    return Promise.all(recordPromises).then((values) => {
      return values;
    });


  }

  parseChromsizes(chrom_names, chrom_sizes) {
    const cumValues = [];
    const chromLengths = {};
    const chrPositions = {};

    let totalLength = 0;

    for (let i = 0; i < chrom_sizes.length; i++) {
      const length = Number(chrom_sizes[i]);
      totalLength += length;

      const newValue = {
        id: i,
        chr: chrom_names[i],
        pos: totalLength - length,
      };

      cumValues.push(newValue);
      chrPositions[newValue.chr] = newValue;
      chromLengths[chrom_names[i]] = length;
    }

    return {
      cumPositions: cumValues,
      chrPositions,
      totalLength,
      chromLengths,
    };
  }

  
 
}

export default SequenceLoader;
