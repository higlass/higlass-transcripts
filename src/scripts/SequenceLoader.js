import { RemoteFile } from "generic-filehandle";

class SequenceLoader {
  //constructor(dataConfig) {
  constructor(fastaUrl, faiUrl) {
    //this.dataConfig = dataConfig;
    //this.trackUid = slugid.nice();

    const { IndexedFasta } = require("@gmod/indexedfasta");

    //this.chromInfo = null;

    // this.chromsizePromise = fetch(dataConfig.chromSizesUrl, {
    //   method: "GET",
    // })
    //   .then((response) => response.text())
    //   .then((chrInfoText) => {
    //     const data = tsvParseRows(chrInfoText);
    //     this.chromInfo = this.parseChromsizesRows(data);
    //   });

    const remoteFA = new RemoteFile(fastaUrl);
    const remoteFAI = new RemoteFile(faiUrl);

    this.sequenceFile = new IndexedFasta({
      fasta: remoteFA,
      fai: remoteFAI,
    });
  }

  // tilesetInfo(callback) {
  //   this.tilesetInfoLoading = true;
  //   return this.chromsizePromise
  //     .then(() => {
  //       this.tilesetInfoLoading = false;

  //       const TILE_SIZE = 1024;
  //       const totalLenth = this.chromInfo.totalLength;
  //       const maxZoom = Math.ceil(
  //         Math.log(totalLenth / TILE_SIZE) / Math.log(2)
  //       );

  //       let retVal = {};

  //       retVal = {
  //         tile_size: TILE_SIZE,
  //         bins_per_dimension: TILE_SIZE,
  //         max_zoom: maxZoom,
  //         max_width: TILE_SIZE * 2 ** maxZoom,
  //         min_pos: [0],
  //         max_pos: [totalLenth],
  //       };

  //       if (callback) {
  //         callback(retVal);
  //       }

  //       return retVal;
  //     })
  //     .catch((err) => {
  //       this.tilesetInfoLoading = false;

  //       if (callback) {
  //         callback({
  //           error: `Error parsing chromsizes: ${err}`,
  //         });
  //       } else {
  //         console.error("Could not fetch tileInfo for sequence track.");
  //       }
  //     });
  // }

  // fetchTilesDebounced(receivedTiles, tileIds) {
  //   const tiles = {};
  //   const zoomLevels = [];
  //   const tilePos = [];
  //   const validTileIds = [];
  //   const tilePromises = [];

  //   for (const tileId of tileIds) {
  //     const parts = tileId.split(".");
  //     const z = parseInt(parts[0], 10);
  //     const x = parseInt(parts[1], 10);

  //     if (Number.isNaN(x) || Number.isNaN(z)) {
  //       console.warn("Invalid tile zoom or position:", z, x);
  //       continue;
  //     }
  //     zoomLevels.push(z);
  //     tilePos.push([x]);
  //     validTileIds.push(tileId);
  //     tilePromises.push(this.tile(z, x));
  //   }

  //   Promise.all(tilePromises).then((values) => {
  //     for (let i = 0; i < values.length; i++) {
  //       const validTileId = validTileIds[i];
  //       tiles[validTileId] = {};
  //       tiles[validTileId].dense = values[i];
  //       tiles[validTileId].zoomLevel = zoomLevels[i];
  //       tiles[validTileId].tilePos = tilePos[i];
  //       tiles[validTileId].tilePositionId = validTileId;
  //     }

  //     receivedTiles(tiles);
  //   });
  //   return tiles;
  // }

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

  // get the sequence for a given tile with optionally an additional number of nuleodides in the beginning
  getTile(z, x, tsInfo, frontExcess = 0) {

    const tileWidth = +tsInfo.max_width / 2 ** +z;

    // get the bounds of the tile
    let minX = tsInfo.min_pos[0] + x * tileWidth - frontExcess;
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
