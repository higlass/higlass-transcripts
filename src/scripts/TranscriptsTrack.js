import { scaleLinear } from "d3-scale";
import slugid from 'slugid';

import {
  initializePixiTexts,
  getContainingExon,
  getTileSequenceOffset,
  exonIntersect,
  getNextExon,
  getAminoAcidsForTile,
  reverseChrCoord,
  getMinusStrandSeq,
} from "./utils";

import SequenceLoader from "./SequenceLoader";

const TranscriptsTrack = (HGC, ...args) => {
  if (!new.target) {
    throw new Error(
      'Uncaught TypeError: Class constructor cannot be invoked without "new"'
    );
  }

  // Services
  const { tileProxy } = HGC.services;

  // Utils
  const { colorToHex, trackUtils, absToChr } = HGC.utils;

  // these are default values that are overwritten by the track's options

  const WHITE_HEX = colorToHex("#ffffff");
  const LIGHT_GREY_HEX = colorToHex("#777777");
  
  const BOXPLOT_DEFAULT_ITEM_RGB_NAME = 'Unknown';
  const BOXPLOT_MIN_SCORE_CUTOFF_FOR_FILL_OPACITY = 300;

  /**
   * Initialize a tile. Pulled out from the track so that it
   * can be modified without having to modify the track
   * object (e.g. in an Observable notebook)
   *
   * @param  {TranscriptsTrack} track   The track object
   * @param  {Object} tile    The tile to render
   * @param  {Object} options The track's options
   */
  function externalInitTile(track, tile, options) {
    const { 
      flipText, 
      fontFamily, 
      labelFontSize, 
      labelFontWeight, 
      highlightTranscriptType,
      highlightTranscriptLabelFontWeight, 
      maxTexts,
    } = options;

    // create texts
    tile.texts = {};
    tile.textWidths = {};
    tile.textHeights = {};

    tile.rectGraphics = new HGC.libraries.PIXI.Graphics();
    tile.rectMaskGraphics = new HGC.libraries.PIXI.Graphics();
    tile.codonSeparatorGraphics = new HGC.libraries.PIXI.Graphics();
    tile.codonTextGraphics = new HGC.libraries.PIXI.Graphics();
    tile.labelBgGraphics = new HGC.libraries.PIXI.Graphics();
    tile.labelGraphics = new HGC.libraries.PIXI.Graphics();

    tile.graphics.addChild(tile.rectGraphics);
    tile.graphics.addChild(tile.rectMaskGraphics);
    tile.graphics.addChild(tile.codonSeparatorGraphics);
    tile.graphics.addChild(tile.codonTextGraphics);
    tile.graphics.addChild(tile.labelBgGraphics);
    tile.graphics.addChild(tile.labelGraphics);

    tile.rectGraphics.mask = tile.rectMaskGraphics;

    if (!tile.tileData.sort) return;

    tile.tileData.sort((a, b) => b.importance - a.importance);

    tile.tileData.forEach((td, i) => {
      const ts = td.fields;
      const tsFormatted = track.formatTranscriptData(ts);

      if (!tsFormatted) return;

      const transcriptName = tsFormatted.transcriptName;
      const transcriptId = tsFormatted.transcriptId;
      const strand = tsFormatted.strand;
      const isLongestIsoform = tsFormatted.isLongestIsoform;
      const isApprisPrincipalIsoform = tsFormatted.isApprisPrincipalIsoform;

      td["transcriptId"] = transcriptId;
      td["transcriptName"] = transcriptName;

      // don't draw texts for the latter entries in the tile
      if (i >= maxTexts) return;

      const text = new HGC.libraries.PIXI.Text(transcriptName, {
        fontSize: `${labelFontSize}px`,
        fontWeight: `${labelFontWeight}`,
        fontFamily,
        fill: track.colors["labelFont"],
      });

      text.interactive = true;

      if (flipText) text.scale.x = -1;

      text.anchor.x = 0;
      text.anchor.y = 0.5;
      text.visible = false;

      switch (highlightTranscriptType) {
        case "none":
          break;
        case "longestIsoform":
          text.style.fontWeight = (isLongestIsoform) ? highlightTranscriptLabelFontWeight : labelFontWeight;
          break;
        case "apprisPrincipalIsoform":
          text.style.fontWeight = (isApprisPrincipalIsoform) ? highlightTranscriptLabelFontWeight : labelFontWeight;
          break;
        default:
          throw new Error(
            'Uncaught TypeError: Unknown highlightTranscriptType option (transcript label text)'
          );
      }

      tile.texts[transcriptId] = text; // index by transcriptName
      tile.texts[transcriptId].strand = strand;
      tile.labelGraphics.addChild(text);
    });

    loadAminoAcidData(track, tile);

    tile.initialized = true;
  }

  function loadAminoAcidData(track, tile) {
    if (
      track.zoomLevel !== track.tilesetInfo.max_zoom ||
      track.sequenceLoader === undefined
    ) {
      return;
    }

    tile.aaInfo = {};
    tile.aaInfo["exonOffsets"] = {};
    tile.aaInfo["nucSequences"] = {};
    tile.aaInfo["aminoAcids"] = {};
    tile.aaInfo["tileOffset"] = 0;

    const chromSizes = track.tilesetInfo.chrom_sizes.split("\t").map((x) => +x);
    const chromNames = track.tilesetInfo.chrom_names.split("\t");

    const chromInfo = track.sequenceLoader.parseChromsizes(
      chromNames,
      chromSizes
    );

    const tileId = +tile.tileId.split(".")[1];
    const zoomLevel = +tile.tileId.split(".")[0];

    const tileSequence = track.sequenceLoader.getTile(
      zoomLevel,
      tileId,
      track.tilesetInfo
    );

    // get the bounds of the tile
    const tileWidth = +track.tilesetInfo.max_width / 2 ** zoomLevel;
    const minX = track.tilesetInfo.min_pos[0] + tileId * tileWidth; // abs coordinates
    const maxX = track.tilesetInfo.min_pos[0] + (tileId + 1) * tileWidth;

    const minXlocOrig = +absToChr(minX, chromInfo)[1];
    const maxXlocOrig = +absToChr(maxX, chromInfo)[1];

    // console.log("Tile bounds abs", minX, maxX);
    // console.log("Tile bounds chr", minXlocOrig, maxXlocOrig);
    // console.log(chromInfo)

    // Compute the offsets of each exon, so that we can get codons across exons
    tile.tileData.forEach((td) => {
      const ts = td.fields;
      const transcriptInfo = track.formatTranscriptData(ts);
      const transcriptId = transcriptInfo.transcriptId;
      const chromLength = chromInfo.chromLengths[transcriptInfo.chromName];

      let minXloc = minXlocOrig;
      let maxXloc = maxXlocOrig;

      if (transcriptInfo["startCodonPos"] === "." || transcriptInfo["stopCodonPos"] === ".") return;

      const strand = transcriptInfo["strand"];
      tile.aaInfo["exonOffsets"][transcriptId] = [];
      tile.aaInfo["nucSequences"][transcriptId] = [];
      tile.aaInfo["aminoAcids"][transcriptId] = [];

      // we don't care about the chrOffset here, we can compute the offsets in chr coordinates
      const exonStartsOriginal = transcriptInfo["exonStarts"];
      const exonEndsOriginal = transcriptInfo["exonEnds"];
      let exonStarts = exonStartsOriginal;
      let exonEnds = exonEndsOriginal;
      let startCodonPos = transcriptInfo["startCodonPos"];
      let stopCodonPos = transcriptInfo["stopCodonPos"];

      if (strand === "-") {
        minXloc = reverseChrCoord(maxXloc, chromLength);
        maxXloc = reverseChrCoord(minXloc, chromLength);
        exonStarts = reverseChrCoord(exonEndsOriginal, chromLength);
        exonEnds = reverseChrCoord(exonStartsOriginal, chromLength);
        startCodonPos = reverseChrCoord(startCodonPos, chromLength);
        stopCodonPos = reverseChrCoord(stopCodonPos, chromLength);
      }

      let accumulatedOffset = 0;
      for (let i = 0; i < exonStarts.length; i++) {
        if (exonStarts[i] <= startCodonPos) {
          tile.aaInfo["exonOffsets"][transcriptId].push(0);
        } else {
          const numNucleotidesInPrevExon =
            exonEnds[i - 1] - Math.max(exonStarts[i - 1], startCodonPos);

          const localOffset = (3 - (numNucleotidesInPrevExon % 3)) % 3;
          accumulatedOffset += localOffset;
          const offset = accumulatedOffset % 3;
          tile.aaInfo["exonOffsets"][transcriptId].push(offset);
        }
      }

      tileSequence.then((values) => {
        let seq = values[0];
        if (strand === "-") {
          seq = getMinusStrandSeq(seq);
        }

        const intersection = exonIntersect(
          exonStarts,
          exonEnds,
          startCodonPos,
          stopCodonPos,
          minXloc,
          seq
        );

        // if there are no exons in this tile, stop.
        if (intersection.filter((nuc) => nuc !== ".").length === 0) {
          return;
        }

        tile.aaInfo["nucSequences"][transcriptId].push(intersection);

        let containingExon = null;
        // if the tile starts within an exon, get the sequence offset for the tile
        if (intersection[0] !== ".") {
          containingExon = getContainingExon(exonStarts, exonEnds, minXloc)
            .exon;
          const exonStart = Math.max(startCodonPos, exonStarts[containingExon]);
          const exonOffset =
            tile.aaInfo["exonOffsets"][transcriptId][containingExon];

          tile.aaInfo["tileOffset"] = getTileSequenceOffset(
            exonStart,
            exonOffset,
            minXloc
          );
        } else {
          //Tile stats before start codon -> no tile offset
          if (minXloc < startCodonPos) {
            tile.aaInfo["tileOffset"] = 0;
          } else {
            const nextExon = getNextExon(exonStarts, exonEnds, minXloc);

            tile.aaInfo["tileOffset"] =
              nextExon !== null && minXloc > startCodonPos
                ? tile.aaInfo["exonOffsets"][transcriptId][nextExon.exon]
                : 0;
          }
        }

        getAminoAcidsForTile(
          HGC,
          intersection,
          tile.aaInfo["tileOffset"],
          transcriptInfo.chromName,
          chromLength,
          strand,
          exonStarts,
          exonEnds,
          minXloc,
          track.pixiTexts,
          track.sequenceLoader
        ).then((aa) => {
          tile.aaInfo["aminoAcids"][transcriptId] = aa;
          track.draw();
        });
      });
    });

    return;
  }

  /** Draw the exons within a gene */
  function drawExons(
    track,
    tile,
    transcriptId,
    graphics,
    chrOffset,
    centerY,
    height,
    strand
  ) {
  
    const topY = centerY - height / 2;

    // get the bounds of the tile
    const tileId = +tile.tileId.split(".")[1];
    const zoomLevel = +tile.tileId.split(".")[0]; //track.zoomLevel does not always seem to be up to date
    const tileWidth = +track.tilesetInfo.max_width / 2 ** zoomLevel;
    const tileMinX = track.tilesetInfo.min_pos[0] + tileId * tileWidth; // abs coordinates
    const tileMaxX = track.tilesetInfo.min_pos[0] + (tileId + 1) * tileWidth;
    
    // isoform type
    const transcriptIsLongestIsoform = track.transcriptInfo[transcriptId]["isLongestIsoform"];
    const transcriptIsApprisPrincipalIsoform = track.transcriptInfo[transcriptId]["isApprisPrincipalIsoform"];

    const exonStarts = track.transcriptInfo[transcriptId]["exonStarts"];
    const exonEnds = track.transcriptInfo[transcriptId]["exonEnds"];

    const isProteinCoding =
      track.transcriptInfo[transcriptId]["startCodonPos"] !== "." && track.transcriptInfo[transcriptId]["stopCodonPos"] !== ".";


    const startCodonPos = isProteinCoding
      ? track.transcriptInfo[transcriptId]["startCodonPos"] + chrOffset
      : -1;

    const stopCodonPos = isProteinCoding
      ? track.transcriptInfo[transcriptId]["stopCodonPos"] + chrOffset
      : -1;

    const txStart = track.transcriptInfo[transcriptId]["txStart"] + chrOffset;
    const txEnd = track.transcriptInfo[transcriptId]["txEnd"] + chrOffset;

    let exonOffsetStarts = exonStarts.map((x) => +x + chrOffset);
    let exonOffsetEnds = exonEnds.map((x) => +x + chrOffset);

    // Add start and stop codon to the exon list and distinguish between UTR and coding region later
    if (isProteinCoding) {
      exonOffsetStarts.push(startCodonPos, stopCodonPos);
      exonOffsetEnds.push(startCodonPos, stopCodonPos);

      exonOffsetStarts.sort();
      exonOffsetEnds.sort();
    }

    const xStartPos = track._xScale(txStart + 1);
    const xEndPos = track._xScale(txEnd + 1);

    const width = xEndPos - xStartPos;
    const yMiddle = centerY;

    const polys = [];
    const polysSVG = []; // holds all polygons that need to be drawn for SVG export

    let poly = [];

    // draw track background, if a longest or APPRIS-principal isoform
    switch (track.options.highlightTranscriptType) {
      case "none":
        break;
      case "longestIsoform":
        if (transcriptIsLongestIsoform) {
          graphics.beginFill(track.colors.highlightTrackBackground);
          poly = [
            xStartPos,
            topY,
            xStartPos + width,
            topY,
            xStartPos + width,
            topY + height,
            xStartPos,
            topY + height,
          ];
          graphics.drawPolygon(poly);
          polys.push([xStartPos, xStartPos + width, topY, topY + height]);
          polysSVG.push({
            rect: poly,
            color: track.options.highlightTranscriptTrackBackgroundColor,
            paintOrder: 0
          });
          graphics.endFill();
        }
        break;
      case "apprisPrincipalIsoform":
        if (transcriptIsApprisPrincipalIsoform) {
          graphics.beginFill(track.colors.highlightTrackBackground);
          poly = [
            xStartPos,
            topY,
            xStartPos + width,
            topY,
            xStartPos + width,
            topY + height,
            xStartPos,
            topY + height,
          ];
          graphics.drawPolygon(poly);
          polys.push([xStartPos, xStartPos + width, topY, topY + height]);
          polysSVG.push({
            rect: poly,
            color: track.options.highlightTranscriptTrackBackgroundColor,
            paintOrder: 0
          });
          graphics.endFill();
        }
        break;  
      default:
        throw new Error(
          'Uncaught TypeError: Unknown highlightTranscriptType option (track background fill)'
        );
    }

    const itemRgbIndex = track.transcriptInfo[transcriptId]["itemRgbIndex"];
    const itemRgb = track.transcriptInfo[transcriptId]["itemRgb"];
    const itemRgbColorMode = (itemRgb === "protein_coding") ? "default" : "custom";
    const itemRgbForSVG = `rgb(${itemRgb})`;
    const itemRgbTriplet = itemRgb.split(',');
    const itemRgbFill = HGC.libraries.PIXI.utils.rgb2hex([
      itemRgbTriplet[0] / 255.0,
      itemRgbTriplet[1] / 255.0,
      itemRgbTriplet[2] / 255.0
    ]);

    // draw the middle line
    graphics.beginFill((itemRgbIndex !== -1) ? itemRgbFill : track.colors.intron);
    switch (track.options.blockStyle) {
      case "directional":
        poly = [
          xStartPos,
          yMiddle - 1,
          xStartPos + width,
          yMiddle - 1,
          xStartPos + width,
          yMiddle + 1,
          xStartPos,
          yMiddle + 1,
        ];
        break;
      case "UCSC-like": 
        poly = [
          xStartPos,
          yMiddle - 0.5,
          xStartPos + width,
          yMiddle - 0.5,
          xStartPos + width,
          yMiddle + 0.5,
          xStartPos,
          yMiddle + 0.5,
        ];
        break;
      default:
        throw new Error(
          'Uncaught TypeError: Unknown blockStyle option (middle line)'
        );
    }
    graphics.drawPolygon(poly);
    graphics.endFill();

    // For mouseOver
    polys.push([xStartPos, xStartPos + width, topY, topY + height]);

    // For SVG export
    polysSVG.push({
      rect: poly,
      color: (itemRgbIndex !== -1) ? itemRgbForSVG : track.colors.intronHEX,
      colorMode: itemRgbColorMode,
      paintOrder: 0
    });

    // draw the actual exons
    for (let j = 0; j < exonOffsetStarts.length; j++) {      
      const exonStart = exonOffsetStarts[j];
      const exonEnd = exonOffsetEnds[j];

      // if the exon has no overlap with the tile, go on
      if (exonEnd < tileMinX || exonStart > tileMaxX) {
        continue;
      }

      const isNonCodingOrUtr =
        !isProteinCoding ||
        (strand === "+" &&
          (exonEnd <= startCodonPos || exonStart >= stopCodonPos)) ||
        (strand === "-" &&
          (exonStart >= startCodonPos || exonEnd <= stopCodonPos));

      const colorUsed = (itemRgbIndex !== -1) ? itemRgbFill : isNonCodingOrUtr ? track.colors.utr : track.colors[strand];
      const colorUsedSVG = (itemRgbIndex !== -1) ? itemRgbForSVG : isNonCodingOrUtr ? track.options.utrColor : track.colors[strand+"HEX"];

      // graphics.beginFill(colorUsed);
      const xStart = track._xScale(exonStart + 1);
      const localWidth = Math.max(
        2,
        track._xScale(exonEnd + 1) - track._xScale(exonStart + 1)
      );

      let minX = xStartPos;
      let maxX = xEndPos;
      let localPoly = null;
      let localRect = null; // without direction for mouseOver
      let chevronPoly = null;
      const chevronYPad = height / 4;
      const chevronWidth = 6;
      const chevronWidthOffset = chevronWidth / 2;
      const chevronThreshold = chevronWidth * 6;

      // for use with "UCSC-like" blockStyle type
      let blockRect = [];
      let blockPoly = [];

      const utrYOffset = (isNonCodingOrUtr) ? height / 4 : 0;
      const utrXOffset = (isNonCodingOrUtr) ? 1 : 0;

      if (strand === "+") {
        const rectStartX = Math.min(xStart, maxX);
        const rectStartX2 = Math.max(rectStartX - 5, xStartPos);
        const rectEndX = Math.min(xStart + localWidth, maxX);
        const rectEndX2 = Math.max(rectEndX - 5, xStartPos);

        localPoly = [
          rectStartX,  topY,
          rectEndX2,   topY,
          rectEndX,    topY + height / 2,
          rectEndX2,   topY + height,
          rectStartX2, topY + height,
          rectStartX,  topY + height / 2,
          rectStartX2, topY,
        ];

        localRect = [rectStartX, rectEndX, topY, topY + height];

        blockRect = [
          rectStartX + utrXOffset, topY + utrYOffset, 
          localWidth - utrXOffset, height - 2 * utrYOffset
        ];
        blockPoly = [
          rectStartX + utrXOffset, topY + utrYOffset, 
          rectStartX + utrXOffset, topY + height - utrYOffset, 
          rectEndX - utrXOffset,   topY + height - utrYOffset, 
          rectEndX - utrXOffset,   topY + utrYOffset,
          rectStartX + utrXOffset, topY + utrYOffset, 
        ];
        
        if (j < exonOffsetStarts.length - 1) {
          const forIntronStart = exonEnd + 1;
          const forIntronEnd = exonOffsetStarts[j + 1] - 1;
          const forIntronXStart = track._xScale(forIntronStart + 1);
          const forIntronLocalWidth = Math.max(
            0,
            track._xScale(forIntronEnd + 1) - track._xScale(forIntronStart + 1)
          );
          if (forIntronLocalWidth > chevronThreshold) {
            const forChevronXStart = forIntronXStart + (forIntronLocalWidth / 2) - chevronWidthOffset;
            chevronPoly = [
              forChevronXStart,                yMiddle - chevronYPad,
              forChevronXStart + chevronWidth, yMiddle,
              forChevronXStart,                yMiddle + chevronYPad,
            ];
            graphics.beginFill((itemRgbIndex !== -1) ? itemRgbFill : track.colors.intron);
            graphics.drawPolygon(chevronPoly);
            graphics.endFill();
            polysSVG.push({
              rect: chevronPoly,
              color: (itemRgbIndex !== -1) ? colorUsedSVG : track.colors.intronHEX,
              colorMode: itemRgbColorMode,
              paintOrder: 0
            });
          }
        }
      } 
      else {
        const rectStartX = Math.max(xStart, minX);
        const rectStartX2 = Math.min(rectStartX + 5, xEndPos);
        const rectEndX = Math.min(Math.max(xStart + localWidth, minX),xEndPos);
        const rectEndX2 = Math.min(rectEndX + 5, xEndPos);

        localPoly = [
          rectStartX,  topY + height / 2,
          rectStartX2, topY,
          rectEndX2,   topY,
          rectEndX,    topY + height / 2,
          rectEndX2,   topY + height,
          rectStartX2, topY + height,
          rectStartX,  topY + height / 2,
        ];

        localRect = [rectStartX, rectEndX, topY, topY + height];

        blockRect = [
          rectStartX + utrXOffset, topY + utrYOffset, 
          localWidth - utrXOffset, height - 2 * utrYOffset
        ];
        blockPoly = [
          rectStartX + utrXOffset, topY + utrYOffset, 
          rectStartX + utrXOffset, topY + height - utrYOffset, 
          rectEndX - utrXOffset,   topY + height - utrYOffset, 
          rectEndX - utrXOffset,   topY + utrYOffset,
          rectStartX + utrXOffset, topY + utrYOffset,
        ];

        if (j < exonOffsetStarts.length - 1) {
          const revIntronStart = exonEnd + 1;
          const revIntronEnd = exonOffsetStarts[j + 1] - 1;
          const revIntronXStart = track._xScale(revIntronStart + 1);
          const revIntronLocalWidth = Math.max(
            0,
            track._xScale(revIntronEnd + 1) - track._xScale(revIntronStart + 1)
          );
          if (revIntronLocalWidth > chevronThreshold) {
            const revChevronXStart = revIntronXStart + (revIntronLocalWidth / 2) + chevronWidthOffset;
            chevronPoly = [
              revChevronXStart,                yMiddle - chevronYPad,
              revChevronXStart - chevronWidth, yMiddle,
              revChevronXStart,                yMiddle + chevronYPad,
            ];
            graphics.beginFill((itemRgbIndex !== -1) ? itemRgbFill : track.colors.intron);
            graphics.drawPolygon(chevronPoly);
            graphics.endFill();
            polysSVG.push({
              rect: chevronPoly,
              color: (itemRgbIndex !== -1) ? colorUsedSVG : track.colors.intronHEX,
              colorMode: itemRgbColorMode,
              paintOrder: 0
            });
          }
        }
      }

      graphics.beginFill(colorUsed);
      switch (track.options.blockStyle) {

        case "directional":
          graphics.drawPolygon(localPoly);
          graphics.endFill();
          polys.push(localRect);
          // For SVG export
          polysSVG.push({
            rect: localPoly,
            color: colorUsedSVG,
            colorMode: itemRgbColorMode,
            paintOrder: 1
          });
          break;

        case "UCSC-like":
          graphics.drawRect(...blockRect);
          graphics.endFill();
          polys.push(localRect);
          // For SVG export
          polysSVG.push({
            rect: blockPoly,
            color: colorUsedSVG,
            colorMode: itemRgbColorMode,
            paintOrder: 1
          });
          break;

        default:
          throw new Error(
            'Uncaught TypeError: Unknown blockStyle option (exon)'
          );
      }      
    }

    const polysForMouseOver = polys.map((x) => [
      x,
      track.transcriptInfo[transcriptId],
    ]);

    tile.allExonsForMouseOver = tile.allExonsForMouseOver.concat(
      polysForMouseOver
    );

    tile.allExonsForSVG = tile.allExonsForSVG.concat(
      polysSVG
    );

    return;
  }

  function renderTranscriptExons(
    transcripts,
    track,
    tile,
    centerY,
    height,
    strandSpacing
  ) {
    transcripts.forEach((transcript) => {
      const transcriptInfo = transcript.fields;
      const chrOffset = +transcript.chrOffset;

      const transcriptId = track.transcriptId(transcriptInfo);

      if (!transcript) return;

      if (!transcriptId) return;

      if (!track.transcriptInfo) return;

      if (!track.transcriptInfo[transcriptId]) return;

      if (track.areTranscriptsHidden && track.transcriptInfo[transcriptId].displayOrder !== 0){
        return;
      };

      if(!track.transcriptInfo[transcriptId]){
        return;
      }

      let centerYOffset =
        track.transcriptInfo[transcriptId].displayOrder *
        (height + strandSpacing);

      if (track.options.showToggleTranscriptsButton) {
        centerYOffset += track.toggleButtonHeight;
      }

      switch (track.options.blockStyle) {
        case "directional":
        case "UCSC-like":
          drawExons(
            track,
            tile,
            transcriptId,
            tile.rectGraphics, //graphics,
            chrOffset,
            centerY + centerYOffset,
            height,
            transcript.strand || transcript.fields[5]
          );
          break;
        case "boxplot":
          drawBoxplotElement(
            track,
            tile,
            transcriptId,
            tile.rectGraphics, //graphics,
            chrOffset,
            centerY + centerYOffset,
            height
          );
          break;
        default:
          throw new Error(
            'Uncaught TypeError: Unknown blockStyle option (drawExon/drawBoxplot)'
          );
      }
    });
  }

  function drawBoxplotElement(
    track,
    tile,
    transcriptId,
    graphics,
    chrOffset,
    centerY,
    height
  ) {
    const topY = centerY - height / 2;

    const txStart = track.transcriptInfo[transcriptId]["txStart"] + chrOffset;
    const txEnd = track.transcriptInfo[transcriptId]["txEnd"] + chrOffset;

    const xStartPos = track._xScale(txStart + 1);
    const xEndPos = track._xScale(txEnd + 1);

    let width = xEndPos - xStartPos;
    const yMiddle = centerY;

    const polys = [];
    const polysSVG = []; // holds all polygons that need to be drawn for SVG export

    let poly = [];

    const boxplotElementFillTriplet = track.transcriptInfo[transcriptId]["itemRgb"].split(',');
    const boxplotElementFill = HGC.libraries.PIXI.utils.rgb2hex([
      boxplotElementFillTriplet[0] / 255.0,
      boxplotElementFillTriplet[1] / 255.0,
      boxplotElementFillTriplet[2] / 255.0
    ]);

    const boxplotElementScore = +track.transcriptInfo[transcriptId]["importance"];
    let boxplotElementFillOpacity =
      (boxplotElementScore < BOXPLOT_MIN_SCORE_CUTOFF_FOR_FILL_OPACITY
        ? BOXPLOT_MIN_SCORE_CUTOFF_FOR_FILL_OPACITY
        : boxplotElementScore) / 1000.0;

    if (Number.isNaN(boxplotElementFillOpacity)) { 
      console.warn(`Element lacks importance/score data`);
      boxplotElementFillOpacity = 1.0;
    };
    const boxplotElementFillAsRgb = `rgb(${boxplotElementFillTriplet[0]},${boxplotElementFillTriplet[1]},${boxplotElementFillTriplet[2]},${boxplotElementFillOpacity})`;

    // begin drawing
    graphics.beginFill(boxplotElementFill, boxplotElementFillOpacity);
    if (width > 2) {
      // draw the middle line
      poly = [
        xStartPos,
        yMiddle - 1,
        xStartPos + width,
        yMiddle - 1,
        xStartPos + width,
        yMiddle + 1,
        xStartPos,
        yMiddle + 1,
      ];
      graphics.drawPolygon(poly);
      polysSVG.push({
        rect: poly,
        color: boxplotElementFillAsRgb,
        paintOrder: 0
      });
  
      // draw start cap
      poly = [
        xStartPos,
        topY,
        xStartPos + 1,
        topY,
        xStartPos + 1,
        topY + height,
        xStartPos,
        topY + height,
      ];
      graphics.drawPolygon(poly);
      polysSVG.push({
        rect: poly,
        color: boxplotElementFillAsRgb,
        paintOrder: 0
      });

      // end cap
      poly = [
        xEndPos,
        topY,
        xEndPos - 1,
        topY,
        xEndPos - 1,
        topY + height,
        xEndPos,
        topY + height,
      ];
      graphics.drawPolygon(poly);
      polysSVG.push({
        rect: poly,
        color: boxplotElementFillAsRgb,
        paintOrder: 0
      });
      
      // draw box
      for (let j = 1; j < track.transcriptInfo[transcriptId]["boxBlockCount"] - 1; j++) { 
        const boxStart = txStart + track.transcriptInfo[transcriptId]["boxBlockStarts"][j];
        const boxEnd = boxStart + track.transcriptInfo[transcriptId]["boxBlockLengths"][j];
        const xBoxStartPos = track._xScale(boxStart + 1);
        const xBoxEndPos = track._xScale(boxEnd + 1);
        poly = [
          xBoxStartPos,
          topY + height / 4,
          xBoxEndPos,
          topY + height / 4,
          xBoxEndPos,
          topY + 3 * height / 4,
          xBoxStartPos,
          topY + 3 * height / 4,
        ];
        graphics.drawPolygon(poly);
        polysSVG.push({
          rect: poly,
          color: boxplotElementFillAsRgb,
          paintOrder: 0
        });
      }
  
      // draw midpoint
      const startCodonPos = +track.transcriptInfo[transcriptId]["startCodonPos"] + chrOffset;
      const stopCodonPos = +track.transcriptInfo[transcriptId]["stopCodonPos"] + chrOffset;
      const xMidpointStartPos = track._xScale(startCodonPos - 1);
      const xMidpointEndPos = track._xScale(stopCodonPos + 1);
      poly = [
        xMidpointStartPos,
        topY,
        xMidpointEndPos,
        topY,
        xMidpointEndPos,
        topY + height,
        xMidpointStartPos,
        topY + height,
      ];
      graphics.drawPolygon(poly);
      polysSVG.push({
        rect: poly,
        color: boxplotElementFillAsRgb,
        paintOrder: 0
      });
    }
    else {
      // draw vertical line only
      width = 3;
      poly = [
        xStartPos,
        topY,
        xStartPos + width,
        topY,
        xStartPos + width,
        topY + height,
        xStartPos,
        topY + height,
      ];
      graphics.drawPolygon(poly);
      polysSVG.push({
        rect: poly,
        color: boxplotElementFillAsRgb,
        paintOrder: 0
      });
    }
    graphics.endFill();


    // For mouseOver
    polys.push([xStartPos, xStartPos + width, topY, topY + height]);

    const polysForMouseOver = polys.map((x) => [
      x,
      track.transcriptInfo[transcriptId],
    ]);
    
    tile.allExonsForMouseOver = tile.allExonsForMouseOver.concat(
      polysForMouseOver
    );

    tile.allExonsForSVG = tile.allExonsForSVG.concat(
      polysSVG
    );

    return;
  }

  /** Create a preventing this track from drawing outside of its
   * visible area
   */
  function renderMask(track, tile) {
    const { tileX, tileWidth } = trackUtils.getTilePosAndDimensions(
      track.tilesetInfo,
      tile.tileId
    );

    tile.rectMaskGraphics.clear();

    const randomColor = Math.floor(Math.random() * 16 ** 6);
    tile.rectMaskGraphics.beginFill(randomColor, 0.3);

    const x = track._xScale(tileX);
    const y = 0;
    const width = (track.isVisible) ? track._xScale(tileX + tileWidth) - track._xScale(tileX) : 0;
    const height = (track.isVisible) ? track.dimensions[1] : 0;
    tile.rectMaskGraphics.drawRect(x, y, width, height);
  }

  const toggleBtnHover = (event, track, overOrOut) => {
    if (overOrOut === "over") {
      track.buttons["pToggleButton"].children[0].alpha = 0.8;
      document.body.style.cursor = "pointer"; // I guess that's not very elegant
    } else if (overOrOut === "out") {
      track.buttons["pToggleButton"].children[0].alpha = 0.5;
      document.body.style.cursor = "default";
    }
    requestAnimationFrame(track.animate);
  };

  const toggleBtnClick = (event, track) => {
    if (!track.areTranscriptsHidden) {
      const numNonVisibleTranscripts = Object.keys(track.transcriptInfo).length - track.transcriptPositionInfo[0].length;
      const numTranscripts = Math.max(0, numNonVisibleTranscripts) ;
      track.buttons["pToggleButton"].children[0].text = "SHOW " + numTranscripts + " MORE TRANSCRIPTS...";
      track.areTranscriptsHidden = true;
    } else {
      track.buttons["pToggleButton"].children[0].text = "SHOW FEWER TRANSCRIPTS...";
      track.areTranscriptsHidden = false;
    }

    track.pubSub.publish("trackDimensionsModified", {
      height: track.computeTrackHeight(),
      resizeParentDiv: true,
      trackId: track.trackId,
      viewId: track.viewId,
    });
  };

  function renderToggleBtn(track) {
    if (
      !track.options.showToggleTranscriptsButton ||
      track.hasToggleBtnBeenRendered
    ) {
      return;
    }
    const pToggleButton = track.buttons["pToggleButton"];
    track.pForeground.removeChildren();
    track.pForeground.addChild(pToggleButton);
    pToggleButton.clear();
    pToggleButton.removeChildren();

    const text = new HGC.libraries.PIXI.Text("SHOW FEWER TRANSCRIPTS...", {
      fontSize: `10px`,
      fontFamily: track.options.fontFamily,
      fontWeight: "500",
      fill: track.colors["black"],
    });
    text.interactive = true;
    text.buttonMode = true;

    text.mouseover = (evt) => toggleBtnHover(evt, track, "over");
    text.mouseout = (evt) => toggleBtnHover(evt, track, "out");
    text.pointerup = (evt) => toggleBtnClick(evt, track);

    text.alpha = 0.5;
    text.anchor.x = 0.5;
    text.anchor.y = 0.5;
    text.position.x = track.dimensions[0] / 2;
    text.position.y = track.toggleButtonHeight / 2;

    pToggleButton.addChild(text);
    track.hasToggleBtnBeenRendered = true;
  }

  class TranscriptsTrackClass extends HGC.tracks
    .HorizontalGeneAnnotationsTrack {
    constructor(context, options) {
      super(context, options);
      const { animate, viewUid } = context;

      this.trackId = this.id;
      this.viewId = viewUid;

      this.animate = animate;

      this.options = options;
      this.initOptions();

      this.numTranscriptRows = 0;

      this.trackHeight = 0;
      this.trackHeightOld = 0;

      this.areTranscriptsHidden = this.options.startCollapsed;
      this.areCodonsShown = false;

      this.transcriptInfo = {};
      this.transcriptSequences = {};

      // Initialize various button infos
      this.buttons = {};

      const pToggleButton = new HGC.libraries.PIXI.Graphics();
      pToggleButton.interactive = true;
      pToggleButton.buttonMode = true;
      this.buttons["pToggleButton"] = pToggleButton;
    }

    initOptions() {
      this.fontSize = +this.options.fontSize;
      this.transcriptHeight = +this.options.transcriptHeight;

      this.transcriptSpacing = +this.options.transcriptSpacing;
      this.geneStrandHSpacing = this.transcriptSpacing / 2;
      this.transcriptHHeight = this.transcriptHeight / 2;

      this.trackHeightAdjustment = this.options.trackHeightAdjustment === "automatic";

      this.toggleButtonHeight = 26;
      // controls when the abbreviated codon text are displayed
      this.minCodonDistance = 15;

      this.codonTextOptions = {
        fontSize: `${this.fontSize * 2}px`,
        fontFamily: this.options.fontFamily,
        fill: WHITE_HEX,
        fontWeight: "bold",
      };

      if (typeof this.options.sequenceData !== "undefined") {
        this.sequenceLoader = new SequenceLoader(
          this.options.sequenceData.fastaUrl,
          this.options.sequenceData.faiUrl
        );
        if (!this.pixiTexts) {
          this.pixiTexts = initializePixiTexts(this.codonTextOptions, HGC);
        }
      }

      this.trackMargin = {
        top: 0,
        left: 0,
        bottom: 10,
        right: 0,
      };
      if (typeof this.options.trackMargin !== "undefined") {
        this.trackMargin = this.options.trackMargin;
      }

      this.isVisible = this.options.isVisible;

      this.colors = {};
      this.colors["+"] = colorToHex(this.options.plusStrandColor);
      this.colors["-"] = colorToHex(this.options.minusStrandColor);
      this.colors["+HEX"] = this.options.plusStrandColor;
      this.colors["-HEX"] = this.options.minusStrandColor;
      this.colors["utr"] = colorToHex(this.options.utrColor);
      this.colors["labelFont"] = colorToHex(this.options.labelFontColor);
      this.colors["black"] = colorToHex("#000000");
      this.colors["intron"] = colorToHex("#CFCFCF");
      this.colors["intronHEX"] = "#CFCFCF";
      this.colors["labelBackgroundPlus"] = colorToHex(this.options.labelBackgroundPlusStrandColor);
      this.colors["labelBackgroundMinus"] = colorToHex(this.options.labelBackgroundMinusStrandColor);
      this.colors["labelStrokePlus"] = colorToHex(this.options.labelStrokePlusStrandColor);
      this.colors["labelStrokeMinus"] = colorToHex(this.options.labelStrokeMinusStrandColor);
      this.colors["background"] = colorToHex(this.options.backgroundColor);
      this.colors["highlightLabelBackground"] = colorToHex(this.options.highlightTranscriptLabelBackgroundColor);
      this.colors["highlightTrackBackground"] = colorToHex(this.options.highlightTranscriptTrackBackgroundColor);
    }

    initTile(tile) {
      externalInitTile(this, tile, {
        flipText: this.flipText,
        fontFamily: this.options.fontFamily,
        labelFontSize: this.options.labelFontSize,
        labelFontWeight: this.options.labelFontWeight,
        highlightTranscriptType: this.options.highlightTranscriptType,
        highlightTranscriptLabelFontWeight: this.options.highlightTranscriptLabelFontWeight,
        maxTexts: this.options.maxTexts,
        blockStyle: this.options.blockStyle,
      });

      // We have to rerender everything since the vertical position
      // of the tracks might have changed accross tiles
      (this.options) && this.rerender(this.options, true);
    }

    /** cleanup */
    destroyTile(tile) {
      tile.rectGraphics.destroy();
      tile.rectMaskGraphics.destroy();
      tile.labelGraphics.removeChildren();
      tile.labelGraphics.destroy();
      tile.labelBgGraphics.destroy();
      tile.codonSeparatorGraphics.destroy();
      tile.codonTextGraphics.removeChildren();
      tile.codonTextGraphics.destroy();
      tile.graphics.destroy();
      tile = null;
    }

    computeTrackHeight() {
      if (!this.options.isVisible) return 0;

      let height = 0;

      if (this.areTranscriptsHidden) {
        height = this.toggleButtonHeight + 
          Math.min(1, this.numTranscriptRows) * (this.transcriptHeight + this.transcriptSpacing) + 
          this.trackMargin.top + this.trackMargin.bottom;
      } 
      else {
        const tbh = this.options.showToggleTranscriptsButton
          ? this.toggleButtonHeight
          : 0;

        height =
          this.numTranscriptRows *
            (this.transcriptHeight + this.transcriptSpacing) +
          tbh +
          this.trackMargin.top + this.trackMargin.bottom;
      }

      if ((this.options.minHeight) && (height < +this.options.minHeight)) {
        height = +this.options.minHeight;
      }

      this.trackHeightOld = this.trackHeight;
      this.trackHeight = height;

      return height;
    }

    adjustTrackHeight() {

      this.computeTrackHeight();

      if (!this.options.isVisible) this.trackHeight = 0;

      this.allowResizeParentDiv = (this.options.allowResizeParentDiv) ? this.options.allowResizeParentDiv : true;

      if (this.trackHeightOld === this.trackHeight) {
        return false;
      }

      this.pubSub.publish("trackDimensionsModified", {
        height: this.trackHeight,
        resizeParentDiv: this.allowResizeParentDiv,
        trackId: this.trackId,
        viewId: this.viewId
      });

      if (this.trackHeightOld === 0) this.rerender(this.options, true);

      return true;
    }

    formatTranscriptData(ts) {
      const strand = ts[5];
      let startCodonPos = ts[11] === "." ? "." : (strand === "+" ? +ts[11] - 1 : +ts[11] + 2);
      let stopCodonPos = ts[12] === "." ? "." : (strand === "+" ? +ts[12] + 2 : +ts[12] - 1);
      const exonStarts = ts[9].split(",").map((x) => +x - 1);
      const exonEnds = ts[10].split(",").map((x) => +x);
      const txStart = +ts[1] - 1;
      const txEnd = +ts[2] - 1;
      const strandedStartCodonPos = ((startCodonPos !== ".") && (stopCodonPos !== ".")) ? ((strand === "+") ? ((startCodonPos < stopCodonPos) ? startCodonPos : stopCodonPos) : ((startCodonPos > stopCodonPos) ? startCodonPos : stopCodonPos)) : ".";
      const strandedStopCodonPos = ((startCodonPos !== ".") && (stopCodonPos !== ".")) ? ((strand === "+") ? ((startCodonPos < stopCodonPos) ? stopCodonPos : startCodonPos) : ((startCodonPos > stopCodonPos) ? stopCodonPos : startCodonPos)) : ".";
      const strandedTxStart = (strand === "+") ? ((txStart < txEnd) ? txStart : txEnd) : ((txStart > txEnd) ? txStart : txEnd);
      const strandedTxEnd = (strand === "+") ? ((txStart < txEnd) ? txEnd : txStart) : ((txStart > txEnd) ? txEnd : txStart);
      
      let boxBlockCount = -1;
      let boxBlockLengths = [];
      let boxBlockStarts = [];
      if (this.options.blockStyle === "boxplot") {
        if (!startCodonPos) { startCodonPos = +ts[6] }
        if (!stopCodonPos) { stopCodonPos = +ts[7] }
        boxBlockCount = +ts[9];
        boxBlockLengths = ts[10].split(",").map((x) => +x);
        boxBlockStarts = ts[11].split(",").map((x) => +x);
      }

      function calculateBlocks() {
        let blocks = [];
        const txDiff = parseFloat(txEnd - txStart);
        let blockIdx = 0;
        let exonTypeIdx = 0;
        let intronTypeIdx = 0;
        let fivePrimeUTRTypeIdx = 0;
        let threePrimeUTRTypeIdx = 0;
        switch (strand) {
          case "+": {
            for (let exonIdx = 0; exonIdx < exonStarts.length; exonIdx++) {
              const exonStart = exonStarts[exonIdx];
              const exonEnd = exonEnds[exonIdx] - 1;
              const exonDiff = exonEnd - exonStart;
              const exonStartOffset = parseFloat(exonStart - txStart);
              if ((startCodonPos !== ".") && (stopCodonPos !== ".")) {
                if ((exonIdx === 0) && (exonStart > strandedTxStart) && (exonStart < strandedStartCodonPos) && (exonEnd < strandedStartCodonPos)) {
                  intronTypeIdx += 1;
                  blocks.push({
                    "type" : "Intron",
                    "subtype" : "",
                    "range" : [
                      0.0,
                      exonStartOffset / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : intronTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonEndOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart === strandedTxStart) && (exonStart < strandedStartCodonPos) && (exonEnd < strandedStartCodonPos)) {
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart < strandedStartCodonPos) && (exonEnd >= strandedStartCodonPos)) {
                  const codonDiff = strandedStartCodonPos - exonStart;
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + codonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "",
                    "range" : [
                      (exonStartOffset + codonDiff) / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart > strandedStartCodonPos) && (exonEnd > strandedStartCodonPos) && (exonStart < strandedStopCodonPos) && (exonEnd < strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype": "",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart < strandedStopCodonPos) && (exonEnd >= strandedStopCodonPos)) {
                  const codonDiff = strandedStopCodonPos - exonStart;
                  exonTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + codonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      (exonStartOffset + codonDiff) / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonEnd < strandedTxEnd) && (exonStart > strandedStopCodonPos) && (exonEnd >= strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx < exonStarts.length - 1) {
                    const nextStrandedExonStartOffset = exonStarts[exonIdx + 1] - txStart;
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        nextStrandedExonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonEnd === strandedTxEnd) && (exonStart > strandedStopCodonPos) && (exonEnd >= strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                }
              }
              else {
                // default block type
                exonTypeIdx += 1;
                blocks.push({ 
                  "type" : "Exon", 
                  "subtype" : "",
                  "range" : [ 
                    exonStartOffset / txDiff, 
                    (exonStartOffset + exonDiff) / txDiff 
                  ],
                  "idx" : blockIdx,
                  "typeIdx" : exonTypeIdx,
                  "subtypeIdx" : null,
                });
                blockIdx += 1;
                if (exonIdx < exonStarts.length - 1) {
                  intronTypeIdx += 1;
                  // const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                  // const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - strandedTxStart - 1;
                  const nextStrandedExonEndOffset = exonEnds[exonIdx] - strandedTxStart - 1;
                  blocks.push({
                    "type" : "Intron",
                    "subtype" : "",
                    "range" : [
                      nextStrandedExonEndOffset / txDiff,
                      exonStartOffset / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : intronTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                }
              }
            }
            break;
          }
          case "-": {
            for (let exonIdx = exonStarts.length - 1; exonIdx >= 0; exonIdx--) {
              const exonStart = exonStarts[exonIdx];
              const exonEnd = exonEnds[exonIdx] - 1;
              const exonDiff = exonEnd - exonStart;
              const exonStartOffset = parseFloat(exonStart - txStart);
              if ((startCodonPos !== ".") && (stopCodonPos !== ".")) {
                if ((exonEnd === strandedTxStart) && (exonStart > strandedStartCodonPos) && (exonEnd >= strandedStartCodonPos)) {
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonEnd < strandedTxStart) && (exonStart > strandedStartCodonPos) && (exonEnd >= strandedStartCodonPos)) {
                  if (exonIdx === exonStarts.length - 1) {
                    intronTypeIdx += 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        (exonStartOffset + exonDiff) / txDiff,
                        1.0,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart < strandedStartCodonPos) && (exonEnd >= strandedStartCodonPos)) {
                  exonTypeIdx += 1;
                  fivePrimeUTRTypeIdx += 1;
                  const codonDiff = startCodonPos - exonStart;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "5'UTR",
                    "range" : [
                      (exonStartOffset + codonDiff) / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : fivePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + codonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : exonTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart > strandedStopCodonPos) && (exonEnd > strandedStopCodonPos) && (exonStart < strandedStartCodonPos) && (exonEnd < strandedStartCodonPos)) {
                  exonTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : exonTypeIdx
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart < strandedStopCodonPos) && (exonEnd >= strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  const codonDiff = stopCodonPos - exonStart;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "",
                    "range" : [
                      (exonStartOffset + codonDiff) / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : exonTypeIdx,
                  });
                  blockIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + codonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart > strandedTxEnd) && (exonStart < strandedStopCodonPos) && (exonEnd < strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  if (exonIdx > 0) {
                    intronTypeIdx += 1;
                    const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                    blocks.push({
                      "type" : "Intron",
                      "subtype" : "",
                      "range" : [
                        nextStrandedExonEndOffset / txDiff,
                        exonStartOffset / txDiff,
                      ],
                      "idx" : blockIdx,
                      "typeIdx" : intronTypeIdx,
                      "subtypeIdx" : null,
                    });
                    blockIdx += 1;
                  }
                }
                else if ((exonStart === strandedTxEnd) && (exonStart < strandedStopCodonPos) && (exonEnd < strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                }
                else if ((exonIdx === 0) && (exonStart > strandedTxEnd) && (exonStart < strandedStopCodonPos) && (exonEnd < strandedStopCodonPos)) {
                  exonTypeIdx += 1;
                  threePrimeUTRTypeIdx += 1;
                  blocks.push({
                    "type" : "Exon",
                    "subtype" : "3'UTR",
                    "range" : [
                      exonStartOffset / txDiff,
                      (exonStartOffset + exonDiff) / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : exonTypeIdx,
                    "subtypeIdx" : threePrimeUTRTypeIdx,
                  });
                  blockIdx += 1;
                  intronTypeIdx += 1;
                  blocks.push({
                    "type" : "Intron",
                    "subtype" : "",
                    "range" : [
                      (exonStartOffset + exonDiff) / txDiff,
                      0.0,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : intronTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                }
              }
              else {
                // default block type
                exonTypeIdx += 1;
                blocks.push({ 
                  "type" : "Exon",
                  "subtype" : "",
                  "range" : [ 
                    exonStartOffset / txDiff, 
                    (exonStartOffset + exonDiff) / txDiff 
                  ],
                  "idx" : blockIdx,
                  "typeIdx" : exonTypeIdx,
                });
                blockIdx += 1;
                if (exonIdx > 0) {
                  intronTypeIdx += 1;
                  const nextStrandedExonEndOffset = exonEnds[exonIdx - 1] - txStart - 1;
                  blocks.push({
                    "type" : "Intron",
                    "subtype" : "",
                    "range" : [
                      nextStrandedExonEndOffset / txDiff,
                      exonStartOffset / txDiff,
                    ],
                    "idx" : blockIdx,
                    "typeIdx" : intronTypeIdx,
                    "subtypeIdx" : null,
                  });
                  blockIdx += 1;
                }
              }
            }
            break;
          }
          default: {
            break;
          }
        }
        return blocks;
      }

      const blocks = calculateBlocks();
      const exonCandidates = blocks.filter((d) => ((d.type === "Exon") && (!d.subtype || d.subtype.length === 0))).slice(-1)[0];
      const intronCandidates = blocks.filter((d) => ((d.type === "Intron") && (!d.subtype || d.subtype.length === 0))).slice(-1)[0];
      const fivePrimeUTRCandidates = blocks.filter((d) => (d.subtype === "5'UTR")).slice(-1)[0];
      const threePrimeUTRCandidates = blocks.filter((d) => (d.subtype === "3'UTR")).slice(-1)[0];
      const blockTypeCounts = {
        "exons" : (exonCandidates ? exonCandidates["typeIdx"] : 0),
        "introns" : (intronCandidates ? intronCandidates["typeIdx"] : 0),
        "fivePrimeUTRs" : (fivePrimeUTRCandidates ? fivePrimeUTRCandidates["subtypeIdx"] : 0),
        "threePrimeUTRs" : (threePrimeUTRCandidates ? threePrimeUTRCandidates["subtypeIdx"] : 0),
      };

      let tags = [];
      let isLongestIsoform = false;
      let isApprisPrincipalIsoform = false;
      let itemRgbIndex = -1;

      if (ts.length >= 14) {
        try {
          tags = JSON.parse(ts[13].replace(/\\/gi, '')); // strip quotation mark escapes
          isLongestIsoform = (tags.findIndex(e => e.includes("hg_longest_isoform")) !== -1);
          isApprisPrincipalIsoform = (tags.findIndex(e => e.includes("appris_principal")) !== -1);      
          const itemRgbIndexHit = tags.findIndex(e => (typeof e === 'object') && ('color' in e));
          if (itemRgbIndexHit !== -1) { 
            ts[8] = `${tags[itemRgbIndexHit].color}`; 
          }
          itemRgbIndex = itemRgbIndexHit;
        }
        catch(err) {
          console.warn(`HGC.TranscriptsTrack - ts[13] parsing error ${JSON.stringify(err)}`)
        }
      }

      if (this.options.showHighlightedTranscriptsOnly && (!isLongestIsoform || !isApprisPrincipalIsoform)) {
        return null;
      }

      const rawTranscriptName = `${ts[3]}`;
      let transcriptName = ts[3];
      switch (this.options.blockStyle) {
        case "directional":
        case "UCSC-like":
          break;
        case "boxplot":
          transcriptName = transcriptName.split('|')[0];
          break;
        default:
          throw new Error(
            'Uncaught TypeError: Unknown blockStyle option (transcript name)'
          );
      }

      const result = {
        transcriptId: this.transcriptId(ts),
        transcriptName: transcriptName,
        rawTranscriptName: rawTranscriptName,
        txStart: txStart,
        txEnd: txEnd,
        strand: strand,
        chromName: ts[0],
        biotype: `${ts[8]}`, // https://www.gencodegenes.org/pages/biotypes.html
        itemRgb: (tags.findIndex(e => (typeof e === 'object') && ('color' in e)) !== -1) ? tags[tags.findIndex(e => (typeof e === 'object') && ('color' in e))].color : ts[8], // https://m.ensembl.org/info/website/upload/bed.html
        itemRgbIndex: +tags.findIndex(e => (typeof e === 'object') && ('color' in e)),
        exonStarts: exonStarts,
        exonEnds: exonEnds,
        startCodonPos: startCodonPos,
        stopCodonPos: stopCodonPos,
        importance: +ts[4],
        boxBlockCount: boxBlockCount,
        boxBlockLengths: boxBlockLengths,
        boxBlockStarts: boxBlockStarts,
        blocks: blocks,
        blockTypeCounts: blockTypeCounts,
        tags: tags,
        isLongestIsoform: tags.includes("hg_longest_isoform"),
        isApprisPrincipalIsoform: tags.includes("appris_principal"),
      };

      return result;
    }

    updateTranscriptInfo() {
      // get all visible transcripts
      const visibleTranscriptsObj = {};

      this.visibleAndFetchedTiles().forEach((tile) => {
        tile.tileData.forEach((ts) => {
          visibleTranscriptsObj[ts.transcriptId] = ts.fields;
        });
      });

      const visibleTranscripts = [];
      for (const tsId in visibleTranscriptsObj) {
        if (visibleTranscriptsObj[tsId]) {
          visibleTranscripts.push(visibleTranscriptsObj[tsId]);
        }
      }

      this.transcriptInfo = {};
      this.transcriptPositionInfo = {};

      this.numTranscriptRows = 0;
      visibleTranscripts
        .sort(function (a, b) {
          return +a[1] - b[1];
        })
        .forEach((ts) => {
          const dpo = this.calculateTranscriptRowNumber(
            this.transcriptPositionInfo,
            +ts[1],
            +ts[2]
          );

          if (this.options.maxRows && (dpo >= this.options.maxRows)) {
            return;
          }

          if (this.transcriptPositionInfo[dpo] === undefined) {
            this.transcriptPositionInfo[dpo] = [];
          }

          this.transcriptPositionInfo[dpo].push([+ts[1], +ts[2], ts[3]]);

          const tsFormatted = this.formatTranscriptData(ts); 

          if (!tsFormatted) return;

          const tInfo = {
            transcriptId: tsFormatted.transcriptId,
            transcriptName: tsFormatted.transcriptName,
            rawTranscriptName: tsFormatted.rawTranscriptName,
            txStart: tsFormatted.txStart,
            txEnd: tsFormatted.txEnd,
            strand: tsFormatted.strand,
            chromName: tsFormatted.chromName,
            biotype: tsFormatted.biotype,
            itemRgb: tsFormatted.itemRgb,
            itemRgbIndex: tsFormatted.itemRgbIndex,
            exonStarts: tsFormatted.exonStarts,
            exonEnds: tsFormatted.exonEnds,
            startCodonPos: tsFormatted.startCodonPos,
            stopCodonPos: tsFormatted.stopCodonPos,
            displayOrder: dpo,
            importance: tsFormatted.importance,
            blocks: tsFormatted.blocks,
            blockTypeCounts: tsFormatted.blockTypeCounts,
            tags: tsFormatted.tags,
            isLongestIsoform: tsFormatted.isLongestIsoform,
            isApprisPrincipalIsoform: tsFormatted.isApprisPrincipalIsoform,
            boxBlockCount: tsFormatted.boxBlockCount,
            boxBlockLengths: tsFormatted.boxBlockLengths,
            boxBlockStarts: tsFormatted.boxBlockStarts,
          };
          this.transcriptInfo[tInfo.transcriptId] = tInfo;
        });

      this.numTranscriptRows = Object.keys(this.transcriptPositionInfo).length;

      // Update the button text
      if (this.areTranscriptsHidden && this.buttons["pToggleButton"].children[0] && this.transcriptPositionInfo[0]) {
        const numNonVisibleTranscripts = Object.keys(this.transcriptInfo).length - this.transcriptPositionInfo[0].length;
        const numTranscripts = Math.max(0, numNonVisibleTranscripts) ;
        this.buttons["pToggleButton"].children[0].text = "SHOW " + numTranscripts + " MORE TRANSCRIPTS...";
      } 
    }

    calculateTranscriptRowNumber(transcriptPositionInfo, txStart, txEnd) {
      const numRows = Object.keys(transcriptPositionInfo).length;

      for (let row = 0; row < numRows; row++) {
        let spaceAvailableOnCurrentRow = true;
        transcriptPositionInfo[row].forEach((ts) => {
          const currentTsStart = ts[0];
          const currentTsEnd = ts[1];

          if (
            (currentTsStart <= txStart && txStart <= currentTsEnd) ||
            (currentTsStart <= txEnd && txEnd <= currentTsEnd)
          ) {
            spaceAvailableOnCurrentRow = false;
          }
        });

        if (spaceAvailableOnCurrentRow) {
          return row;
        }
      }

      // If we are here, there are now available space on the existing rows.
      // Add a new one.
      return numRows;
    }

    /*
     * Redraw the track because the options
     * changed
     */
    rerender(options, force) {
      const strOptions = JSON.stringify(options);
      if (!force && strOptions === this.prevOptions) return;

      this.options = options;

      this.initOptions();

      this.prevOptions = strOptions;

      renderToggleBtn(this);

      this.updateTranscriptInfo();

      // Adjusting the track height leads to a full rerender.
      // No need to rerender again

      if (this.trackHeightAdjustment && this.adjustTrackHeight()) return;

      this.visibleAndFetchedTiles().forEach((tile) => {
        this.renderTile(tile);
      });
    }

    drawTile() {}

    transcriptId(transcriptInfo) {
      return `${transcriptInfo[7]}_${transcriptInfo[0]}_${transcriptInfo[1]}_${transcriptInfo[2]}`;
    }

    exonId(transcriptInfo, exonStart, exonEnd) {
      return `${transcriptInfo[7]}_${transcriptInfo[0]}_${transcriptInfo[1]}_${transcriptInfo[2]}_${exonStart}_${exonEnd}`;
    }

    renderTile(tile) {
      if (!tile.initialized) return;

      tile.allExonsForMouseOver = [];
      tile.allExonsForSVG = [];
      // store the scale at while the tile was drawn at so that
      // we only resize it when redrawing
      tile.drawnAtScale = this._xScale.copy();
      tile.rectGraphics.removeChildren();
      tile.rectGraphics.clear();
      tile.labelBgGraphics.clear();

      const plusTranscripts = tile.tileData.filter(
        (td) =>
          td.type !== "filler" && (td.strand === "+" || td.fields[5] === "+")
      );
      const minusTranscripts = tile.tileData.filter(
        (td) =>
          td.type !== "filler" && (td.strand === "-" || td.fields[5] === "-")
      );
      const unstrandedTranscripts = tile.tileData.filter(
        (td) =>
          td.type !== "filler" && (td.strand === "." || td.fields[5] === ".")
      );

      const strandCenterY = this.trackMargin.top + this.transcriptHeight / 2 + this.transcriptSpacing / 2;

      const renderContext = [
        this,
        tile,
        strandCenterY,
        this.transcriptHeight,
        this.transcriptSpacing,
      ];

      // console.log(`there are ${(plusTranscripts.length + minusTranscripts.length + unstrandedTranscripts.length)} transcripts to be rendered in tile ${tile.tileId} | ${tile.remoteId}!`);

      renderTranscriptExons(plusTranscripts, ...renderContext);
      renderTranscriptExons(minusTranscripts, ...renderContext);
      renderTranscriptExons(unstrandedTranscripts, ...renderContext);

      renderMask(this, tile);

      trackUtils.stretchRects(this, [
        (x) => x.rectGraphics,
        (x) => x.rectMaskGraphics,
      ]);
    }

    calculateZoomLevel() {
      // offset by 2 because 1D tiles are more dense than 2D tiles
      // 1024 points per tile vs 256 for 2D tiles

      const xZoomLevel = tileProxy.calculateZoomLevel(
        this._xScale,
        this.tilesetInfo.min_pos[0],
        this.tilesetInfo.max_pos[0]
      );

      let zoomLevel = Math.min(xZoomLevel, this.maxZoom);
      zoomLevel = Math.max(zoomLevel, 0);

      //console.log(zoomLevel, this._xScale.domain())
      return zoomLevel;
    }

    draw() {

      if (!this.isVisible) return;

      trackUtils.stretchRects(this, [
        (x) => x.rectGraphics,
        (x) => x.rectMaskGraphics,
      ]);

      this.areCodonsShown = false;

      this.drawLabels();
      this.drawCodonSeparators();
      this.drawCodonTexts();

      // otherwise codons are not displayed on startup
      requestAnimationFrame(this.animate);

      this.pubSub.publish("trackDimensionsModified", {
        height: this.trackHeight,
        resizeParentDiv: this.allowResizeParentDiv,
        trackId: this.trackId,
        viewId: this.viewId
      });
    }

    drawCodonTexts() {
      Object.values(this.fetchedTiles)
        .filter((tile) => tile.drawnAtScale)
        .forEach((tile) => {
          if (!tile.initialized || !tile.aaInfo || !tile.aaInfo.aminoAcids) {
            return;
          }

          tile.codonTextGraphics.clear();
          tile.codonTextGraphics.removeChildren();
          // if individual codons are too close together, don't draw anything
          let alpha = 1.0;
          const codonWidth = this._xScale(3) - this._xScale(0);
          if (codonWidth < this.minCodonDistance) {
            return;
          } else if (
            codonWidth < this.minCodonDistance + 3 &&
            codonWidth >= this.minCodonDistance
          ) {
            // gracefully fade out
            const alphaScale = scaleLinear()
              .domain([this.minCodonDistance, this.minCodonDistance + 3])
              .range([0, 1])
              .clamp(true);
            alpha = alphaScale(codonWidth);
          }

          const graphics = tile.codonTextGraphics;
          graphics.beginFill(WHITE_HEX);
          graphics.alpha = alpha;

          tile.tileData.forEach((td) => {
            const transcriptId = td.transcriptId;
            if (
              !this.transcriptInfo[transcriptId] ||
              !tile.aaInfo.aminoAcids[transcriptId]
            ){
              return;
            }

            const transcript = this.transcriptInfo[transcriptId];

            if (!transcript) return;

            if (this.areTranscriptsHidden && transcript && transcript.displayOrder !== 0) {
              return;
            }

            const chrOffset = +td.chrOffset;
            const codons = tile.aaInfo.aminoAcids[transcriptId];

            let yMiddle =
              transcript.displayOrder *
                (this.transcriptHeight + this.transcriptSpacing) +
              this.transcriptHeight / 2 +
              this.transcriptSpacing / 2 -
              this.fontSize / 2 -
              1;

            if (this.options.showToggleTranscriptsButton) {
              yMiddle += this.toggleButtonHeight;
            }

            for (var i = 0, n = codons.length; i < n; ++i) {
              const codon = codons[i];

              const startPosX = this._xScale(codon.posStart + chrOffset + 1);
              const endPosX = this._xScale(codon.posEnd + chrOffset + 1);

              // if the codon is not visible, don't bother
              if (endPosX < 0 || startPosX > this.dimensions[0]) {
                continue;
              }

              //const availableSpace = this._xScale((codon.posStart)) - this._xScale((codon.posEnd + 1));
              //console.log(availableSpace);
              if (codonWidth < this.minCodonDistance + 10) {
                const xMiddle =
                  this._xScale(
                    (codon.posStart + codon.posEnd + 1) / 2 + chrOffset + 1
                  ) -
                  codon.widthAbbrev / 2; //(codon.posStart + codon.posEnd + 1) / 2 + chrOffset
                codon.spriteAbbrev.position.x = xMiddle;
                codon.spriteAbbrev.position.y = yMiddle;
                graphics.addChild(codon.spriteAbbrev);
              } else {
                const xMiddle =
                  this._xScale(
                    (codon.posStart + codon.posEnd + 1) / 2 + chrOffset + 1
                  ) -
                  codon.width / 2; //(codon.posStart + codon.posEnd + 1) / 2 + chrOffset
                codon.sprite.position.x = xMiddle;
                codon.sprite.position.y = yMiddle;
                graphics.addChild(codon.sprite);
              }
            }
          });
        });
    }

    drawCodonSeparators() {
      Object.values(this.fetchedTiles)
        .filter((tile) => tile.drawnAtScale)
        .forEach((tile) => {
          if (!tile.initialized || !tile.aaInfo || !tile.aaInfo.aminoAcids)
            return;

          tile.allCodonsForMouseOver = [];
          tile.codonSeparatorGraphics.clear();

          // if individual codons are too close together, don't draw anything
          let alpha = 1.0;
          const codonWidth = this._xScale(3) - this._xScale(0);

          if (codonWidth < this.minCodonDistance) {
            return;
          } else if (
            codonWidth < this.minCodonDistance + 3 &&
            codonWidth >= this.minCodonDistance
          ) {
            // gracefully fade out
            const alphaScale = scaleLinear()
              .domain([this.minCodonDistance, this.minCodonDistance + 3])
              .range([0, 1])
              .clamp(true);
            alpha = alphaScale(codonWidth);
          }

          this.areCodonsShown = true;

          const graphics = tile.codonSeparatorGraphics;
          graphics.beginFill(WHITE_HEX);
          graphics.alpha = alpha;

          tile.tileData.forEach((td) => {
            const transcriptId = td.transcriptId;
            if (
              !this.transcriptInfo[transcriptId] ||
              !tile.aaInfo.aminoAcids[transcriptId]
            ){
              return;
            }
              
            const transcript = this.transcriptInfo[transcriptId];

            if (!transcript) return;

            if (this.areTranscriptsHidden && transcript && transcript.displayOrder !== 0) {
              return;
            }

            const chrOffset = +td.chrOffset;
            const codons = tile.aaInfo.aminoAcids[transcriptId];

            let yMiddle =
              transcript.displayOrder *
                (this.transcriptHeight + this.transcriptSpacing) +
              this.transcriptHeight / 2 +
              this.transcriptSpacing / 2;

            if (this.options.showToggleTranscriptsButton) {
              yMiddle += this.toggleButtonHeight;
            }

            let yMiddleAbove = yMiddle - this.transcriptHeight / 2;
            let yMiddleBelow = yMiddle + this.transcriptHeight / 2;

            for (var i = 0, n = codons.length; i < n; ++i) {
              const codon = codons[i];

              const startPosX = this._xScale(codon.posStart + chrOffset + 1);
              const endPosX = this._xScale(codon.posEnd + chrOffset + 2);

              // if the codon is not visible, don't bother
              if (endPosX < 0 || startPosX > this.dimensions[0]) {
                continue;
              }

              const rectArr = [startPosX, endPosX, yMiddleAbove, yMiddleBelow];
              const codonForMouseOver = [rectArr, codon.aminoAcid];
              tile.allCodonsForMouseOver.push(codonForMouseOver);

              if (
                (transcript.strand === "+" &&
                  codon.posStart === transcript.startCodonPos) ||
                (transcript.strand === "+" &&
                  transcript.exonStarts.includes(codon.posStart)) ||
                (transcript.strand === "-" &&
                  codon.posEnd + 1 === transcript.startCodonPos) ||
                (transcript.strand === "-" &&
                  transcript.exonEnds.includes(codon.posEnd + 1))
              ) {
                continue;
              }

              const caretPos = transcript.strand === "+" ? startPosX : endPosX;
              let rectStartX = caretPos;
              let rectStartX2 = rectStartX;
              let rectEndX = rectStartX;
              let rectEndX2 = rectStartX;

              if (transcript.strand === "+") {
                rectStartX2 = rectStartX2 - 5;
                rectEndX = rectEndX + 2;
                rectEndX2 = rectEndX2 - 3;
              } else {
                rectStartX2 = rectStartX2 + 5;
                rectEndX = rectEndX - 2;
                rectEndX2 = rectEndX2 + 3;
              }

              const localPoly = [
                rectStartX,
                yMiddleAbove,
                rectEndX2,
                yMiddleAbove,
                rectEndX,
                yMiddle,
                rectEndX2,
                yMiddleBelow,
                rectStartX2,
                yMiddleBelow,
                rectStartX,
                yMiddle,
                rectStartX2,
                yMiddleAbove,
              ];

              graphics.drawPolygon(localPoly);
            }
          });
        });
    }

    drawLabels() {
      this.allTexts = [];
      this.allBoxes = [];
      const allTiles = [];

      Object.values(this.fetchedTiles)
        // tile hasn't been drawn properly because we likely got some
        // bogus data from the server
        .filter((tile) => tile.drawnAtScale)
        .forEach((tile) => {
          tile.labelBgGraphics.clear();
          switch (this.options.blockStyle) {
            case "directional":
            case "UCSC-like":
              tile.labelBgGraphics.beginFill(
                typeof this.options.labelBackgroundColor !== "undefined"
                  ? this.colors["labelBackground"]
                  : WHITE_HEX
              );
              break;
            case "boxplot":
              // undefined until we iterate over tile data objects
              break;
            default:
              throw new Error(
                'Uncaught TypeError: Unknown blockStyle option (drawLabels)'
              );
          }

          if (!tile.initialized) return;

          tile.tileData.forEach((td) => {
            // tile probably hasn't been initialized yet
            if (!tile.texts) return;
            
            const transcriptId = td.transcriptId;

            if (this.transcriptInfo[transcriptId] === undefined) return;

            const transcript = this.transcriptInfo[transcriptId];

            if (!transcript) return;
            
            switch (this.options.blockStyle) {
              case "directional":
              case "UCSC-like":
                break;
              case "boxplot":
                const boxplotFillTriplet = transcript.itemRgb.split(',');
                const boxplotFill = PIXI.utils.rgb2hex([
                  boxplotFillTriplet[0] / 255.0,
                  boxplotFillTriplet[1] / 255.0,
                  boxplotFillTriplet[2] / 255.0
                ]);
                tile.labelBgGraphics.beginFill(boxplotFill);                
                break;
              default:
                throw new Error(
                  'Uncaught TypeError: Unknown blockStyle option (drawLabels, within tile data)'
                );
            }

            const text = tile.texts[transcriptId];

            if (!text) return;

            if (!transcript) return;

            if (this.areTranscriptsHidden && transcript && transcript.displayOrder !== 0) {
              text.visible = false;
              return;
            }

            if (!tile.textWidths[transcriptId]) {
              // if we haven't measured the text's width in renderTile, do it now
              // this can occur if the same gene is in more than one tile, so its
              // dimensions are measured for the first tile and not for the second
              const textWidth = text.getBounds().width;
              const textHeight = text.getBounds().height;

              tile.textHeights[transcriptId] = textHeight;
              tile.textWidths[transcriptId] = textWidth;
            }

            const TEXT_MARGIN = 3;
            const CARET_MARGIN = 3;
            const CARET_WIDTH = 5;
            const chrOffset = +td.chrOffset;
            const txStart = transcript["txStart"] + chrOffset;
            const txEnd = transcript["txEnd"] + chrOffset;

            let textYMiddleOffset =
              transcript.displayOrder *
              (this.transcriptHeight + this.transcriptSpacing);

            if (this.options.showToggleTranscriptsButton) {
              textYMiddleOffset += this.toggleButtonHeight;
            }

            let textYMiddle =
              this.transcriptHeight / 2 +
              this.transcriptSpacing / 2 +
              textYMiddleOffset;

            textYMiddle += this.trackMargin.top;

            // take care of label positioning at start or end of transcripts
            if (transcript.strand === "+") {
              text.position.x = Math.max(
                this._xScale(this.xScale().domain()[0]) + TEXT_MARGIN,
                this._xScale(txStart + 1) -
                  tile.textWidths[transcriptId] -
                  2 * TEXT_MARGIN -
                  CARET_MARGIN
              );
            } else {
              text.position.x = Math.max(
                this._xScale(this.xScale().domain()[0]) +
                  TEXT_MARGIN +
                  CARET_WIDTH,
                this._xScale(txStart + 1) -
                  tile.textWidths[transcriptId] -
                  2 * TEXT_MARGIN
              );
            }

            const marginRight =
              transcript.strand === "+"
                ? tile.textWidths[transcriptId] +
                  this.transcriptHeight / 2 +
                  2 * TEXT_MARGIN -
                  CARET_MARGIN
                : tile.textWidths[transcriptId] + TEXT_MARGIN;

            text.position.x = Math.min(
              text.position.x,
              this._xScale(txEnd + 1) - marginRight
            );

            text.position.y = textYMiddle;

            // Determine if the current text should be hidden
            let showText = true;
            if (!transcript) return;
            const dpo = transcript.displayOrder;

            this.transcriptPositionInfo[dpo]
              .filter((ts) => {
                // Check the ones that are left of the current transcript
                return ts[1] < transcript.txStart;
              })
              .forEach((ts) => {
                const endOfTranscript = this._xScale(ts[1] + chrOffset + 1);

                if (endOfTranscript > text.position.x - 4 * TEXT_MARGIN) {
                  showText = false;
                }
                else {
                  showText = true;
                }
              });

            if (showText) {
              text.visible = true;

              let textFill = text.style.fill;
              let boxFillColor = (text.strand === "+") ? this.colors["labelBackgroundPlus"] : this.colors["labelBackgroundMinus"];
              switch (this.options.blockStyle) {
                case "directional":
                case "UCSC-like":
                  switch (this.options.highlightTranscriptType) {
                    case "none":
                      break;
                    case "longestIsoform":
                      boxFillColor = (transcript.isLongestIsoform) ? this.colors.highlightLabelBackground : boxFillColor;
                      break;
                    case "apprisPrincipalIsoform":
                      boxFillColor = (transcript.isApprisPrincipalIsoform) ? this.colors.highlightLabelBackground : boxFillColor;
                      break;
                    default:
                      throw new Error(
                        'Uncaught TypeError: Unknown highlightTranscriptType option (transcript label background)'
                      );
                  }
                  break;
                case "boxplot":
                  // const boxplotFillColorTriplet = transcript.itemRgb.split(',');
                  // const normalizedBoxplotFillColorTriplet = boxplotFillColorTriplet.map((d) => d/255.0);
                  // const boxplotFillColor = PIXI.utils.rgb2hex([
                  //   normalizedBoxplotFillColorTriplet[0], 
                  //   normalizedBoxplotFillColorTriplet[1],
                  //   normalizedBoxplotFillColorTriplet[2]]);
                  // boxFillColor = boxplotFillColor;
                  // // measure luminance to decide label color, ref. https://www.w3.org/TR/WCAG20/ (contrast ratio)
                  // const luminanceThreshold = 0.1; // should be 0.179, but adjusted for component palette
                  // const luminance = function(r, g, b) { return 0.2126 * r + 0.7152 * g + 0.0722 * b; };
                  // textFill = (luminance(...normalizedBoxplotFillColorTriplet) > luminanceThreshold) ? WHITE_HEX : BLACK_HEX;
                  // if (transcript.itemRgb === "255,229,0") { textFill = BLACK_HEX; }
                  textFill = LIGHT_GREY_HEX;
                  break;
                default:
                  throw new Error(
                    'Uncaught TypeError: Unknown blockStyle option (drawLabels, within tile data, B)'
                  );
              }
              text.style.fill = textFill;

              this.allBoxes.push([
                text.position.x - TEXT_MARGIN,
                textYMiddle,
                tile.textWidths[transcriptId] + 2 * TEXT_MARGIN,
                this.transcriptHeight,
                transcript.transcriptName,
                boxFillColor,
              ]);

              this.allTexts.push({
                importance: transcript.importance,
                text,
                caption: transcript.transcriptName,
                strand: transcript.strand,
                isLongestIsoform: transcript.isLongestIsoform,
                isApprisPrincipalIsoform: transcript.isApprisPrincipalIsoform,
              });

              allTiles.push(tile.labelBgGraphics);
            } else {
              text.visible = false;
            }
          });
        });

      switch (this.options.blockStyle) {
        case "directional":
        case "UCSC-like":
          this.renderTextBg(this.allBoxes, this.allTexts, allTiles);
          break;
        case "boxplot":
          //this.renderDirectionlessTextBg(this.allBoxes, this.allTexts, allTiles);
          break;
        default:
          throw new Error(
            'Uncaught TypeError: Unknown blockStyle option (drawLabels, within tile data, B)'
          );
      }
    }

    renderDirectionlessTextBg(allBoxes, allTexts, allTiles) {
      allTexts.forEach((text, i) => {
        if (text.text.visible && allBoxes[i] && allTiles[i]) {
          const [minX, minY, width, height] = allBoxes[i];
          const margin = 1;
          const xl = minX;
          const xlm = xl + margin;
          const xr = minX + width;
          const xrm = xr - margin;
          const yb = minY - height / 2;
          const ybm = yb + margin;
          const yt = minY + height / 2;
          const ytm = yt - margin;

          allTiles[i].beginFill(this.colors["labelStrokePlus"]);

          // direction-less label
          let polyBorder = [xl, yb, xr, yb, xr, yt, xl, yt, xl, yb];

          allTiles[i].drawPolygon(polyBorder);

          allTiles[i].beginFill(allBoxes[i][5]);

          // direction-less label
          let poly = [xlm, ybm, xrm, ybm, xrm, ytm, xlm, ytm, xlm, ybm];
          allTiles[i].drawPolygon(poly);
        }
      });
    }

    renderTextBg(allBoxes, allTexts, allTiles) {
      allTexts.forEach((text, i) => {
        if (text.text.visible && allBoxes[i] && allTiles[i]) {
          const [minX, minY, width, height] = allBoxes[i];
          const margin = 1;
          const xl = minX;
          const xlm = xl + margin;
          const xr = minX + width;
          const xrm = xr - margin;
          const yb = minY - height / 2;
          const ybm = yb + margin;
          const yt = minY + height / 2;
          const ytm = yt - margin;

          if (text.strand === "+") {
            allTiles[i].beginFill(this.colors["labelStrokePlus"]);

            // Directional label
            let polyBorder = [xl, yb, xr, yb, xr + 5, minY, xr, yt, xl, yt];

            allTiles[i].drawPolygon(polyBorder);

            allTiles[i].beginFill(allBoxes[i][5]);

            // Directional label
            let poly = [xlm, ybm, xr, ybm, xrm + 3, minY, xr, ytm, xlm, ytm];
            allTiles[i].drawPolygon(poly);
          } 
          else {
            allTiles[i].beginFill(this.colors["background"]);

            let polyBg = [xl - 5, yb, xr, yb, xr, yt, xl - 5, yt];

            allTiles[i].drawPolygon(polyBg);

            allTiles[i].beginFill(this.colors["labelStrokeMinus"]);
            // Directional label
            let polyBorder = [xl - 5, minY, xl, yb, xr, yb, xr, yt, minX, yt];

            allTiles[i].drawPolygon(polyBorder);

            allTiles[i].beginFill(allBoxes[i][5]);

            // Directional label
            let poly = [
              xlm - 3,
              minY,
              minX,
              ybm,
              xrm,
              ybm,
              xrm,
              ytm,
              minX,
              ytm,
            ];
            allTiles[i].drawPolygon(poly);
          }

        }
      });
    }

    hideTexts(allTexts) {
      allTexts.forEach((text, i) => {
        text.visible = false;
      });
    }

    setPosition(newPosition) {
      super.setPosition(newPosition);
      [this.pMain.position.x, this.pMain.position.y] = this.position;
    }

    setDimensions(newDimensions) {
      this.updateTranscriptInfo();
      // This will rerender all tiles.
      super.setDimensions(newDimensions);
    }

    zoomed(newXScale, newYScale) {
      this.xScale(newXScale);
      this.yScale(newYScale);
      this.refreshTiles();
      this.draw();
    }

    isPointInRectangle(rect, point) {
      if (
        rect[0] < point[0] &&
        rect[1] > point[0] &&
        rect[2] < point[1] &&
        rect[3] > point[1]
      ) {
        return true;
      }
      return false;
    }
    
    getNormPointXWithinLocalRect(pointX, rectLeftX, rectRightX) {
      const rectStartOffset = (rectRightX > rectLeftX) ? rectLeftX : rectRightX;
      const rectDiff = Math.abs(rectRightX - rectLeftX);
      const normPointX = (pointX - rectStartOffset) / rectDiff;
      return (normPointX < 0.0) ? 0.0 : (normPointX > 1.0) ? 1.0 : normPointX;
    }

    formattedBED12HTML(bed12FieldsObj) {
      const chrom = bed12FieldsObj.chrom;
      const start = +bed12FieldsObj.start;
      const end = +bed12FieldsObj.end;
      const id = bed12FieldsObj.id;
      const score = bed12FieldsObj.score;
      const strand = bed12FieldsObj.strand;
      const thickStart = +bed12FieldsObj.thickStart;
      const thickEnd = +bed12FieldsObj.thickEnd;
      const itemRGB = bed12FieldsObj.itemRGB;
      const blockCount = +bed12FieldsObj.blockCount;
      const blockSizes = bed12FieldsObj.blockSizes;
      const blockStarts = bed12FieldsObj.blockStarts;
  
      const idElems = id.split('|');
      const realId = idElems[0];
      const realScorePrecision = 4;
      const realScore =
        idElems.length > 1
          ? Number.parseFloat(idElems[1]).toPrecision(realScorePrecision)
          : score;
  
      const hc = document.getElementsByClassName('higlass-main-content')[0];
      if (hc) {
        hc.style.cursor = 'pointer';
      }
  
      let itemRGBMarkup = '';
      if (this.options.itemRGBMap) {
        const itemRGBName = this.options.itemRGBMap[itemRGB]
          ? this.options.itemRGBMap[itemRGB]
          : BOXPLOT_DEFAULT_ITEM_RGB_NAME;
        itemRGBMarkup = `<div id="bed12-component" style="display:inline-block; position:relative; top:-2px;">
          <svg width="10" height="10">
            <rect width="10" height="10" rx="2" ry="2" style="fill:rgb(${itemRGB});stroke:black;stroke-width:2;" />
          </svg>
          <span style="position:relative; top:1px; font-weight:600;">${itemRGBName}</span>
        </div>`;
      }
  
      /*
        With gratitude to Pierre Lindenbaum @ biostars
      */
      let elementCartoon = '';
      const elementCartoonWidth = 200;
      const elementCartoonGeneHeight = 30;
      const elementCartoonHeight = elementCartoonGeneHeight + 10;
      const elementCartoonMiddle = elementCartoonHeight / 2;
      function pos2pixel(pos) {
        return ((pos - start) / ((end - start) * 1.0)) * elementCartoonWidth;
      }
      if (blockCount > 0) {
        elementCartoon += `<svg width="${elementCartoonWidth}" height="${elementCartoonHeight}">
          <style type="text/css">
            .ticks {stroke:rgb(${itemRGB});stroke-width:1px;fill:none;}
            .gene {stroke:rgb(${itemRGB});stroke-width:1px;fill:none;}
            .translate { fill:rgb(${itemRGB});fill-opacity:1;}
            .exon { fill:rgb(${itemRGB});fill-opacity:1;}
            .score { fill:rgb(${itemRGB});fill-opacity:1;font:bold 12px sans-serif;}
            .id { fill:rgb(${itemRGB});fill-opacity:1;font:bold 12px sans-serif;}
          </style>
          <defs>
            <path id="ft" class="ticks" d="m -3 -3  l 3 3  l -3 3" />
            <path id="rt" class="ticks" d="m 3 -3  l -3 3  l 3 3" />
          </defs>
        `;
        const ecStart = pos2pixel(start);
        const ecEnd = pos2pixel(end);
        elementCartoon += `<line class="gene" x1=${ecStart} x2=${ecEnd} y1=${elementCartoonMiddle} y2=${elementCartoonMiddle} />`;
        const ecThickStart = pos2pixel(thickStart);
        const ecThickEnd = pos2pixel(thickEnd);
        const ecThickY = elementCartoonMiddle - elementCartoonGeneHeight / 4;
        const ecThickHeight = elementCartoonGeneHeight / 2;
        let ecThickWidth = ecThickEnd - ecThickStart;
        if (this.options.isBarPlotLike) {
          ecThickWidth = ecThickWidth !== 1 ? 1 : ecThickWidth;
        }
        let realIdTextAnchor = '';
        if (ecThickStart < 0.15 * elementCartoonWidth) {
          realIdTextAnchor = 'start';
        } else if (
          ecThickStart >= 0.15 * elementCartoonWidth &&
          ecThickStart <= 0.85 * elementCartoonWidth
        ) {
          realIdTextAnchor = 'middle';
        } else {
          realIdTextAnchor = 'end';
        }
        // const realScoreTextAnchor = (ecThickStart < (0.15 * elementCartoonWidth)) ? "start" : (ecThickStart > (0.85 * elementCartoonWidth)) ? "end" : "middle";
        elementCartoon += `<rect class="translate" x=${ecThickStart} y=${ecThickY} width=${ecThickWidth} height=${ecThickHeight} />`;
        const ecLabelDy = '-0.25em';
        elementCartoon += `<text class="id" text-anchor="${realIdTextAnchor}" x=${ecThickStart} y=${ecThickY} dy=${ecLabelDy}>${realId}</text>`;
        if (strand === '+' || strand === '-') {
          const ecStrandHref = strand === '+' ? '#ft' : '#rt';
          for (let i = 0; i < elementCartoonWidth; i += 10) {
            elementCartoon += `<use x=${i} y=${elementCartoonMiddle} href=${ecStrandHref} />`;
          }
        }
        for (let i = 0; i < blockCount; i++) {
          let ecExonStart = pos2pixel(start + +blockStarts[i]);
          const ecExonY = elementCartoonMiddle - elementCartoonGeneHeight / 8;
          let ecExonWidth = pos2pixel(start + +blockSizes[i]);
          const ecExonHeight = elementCartoonGeneHeight / 4;
          if (this.options.isBarPlotLike) {
            if (i === 0) {
              ecExonStart = ecStart;
              ecExonWidth = ecStart + 1;
            } else if (i === blockCount - 1) {
              ecExonStart = ecEnd - 1;
              ecExonWidth = ecEnd;
            }
          }
          elementCartoon += `<rect class="exon" x=${ecExonStart} y=${ecExonY} width=${ecExonWidth} height=${ecExonHeight} />`;
        }
  
        elementCartoon += '</svg>';
      }
  
      let intervalMarkup = `${chrom}:${start}-${end}`;
      if (strand === '+' || strand === '-') {
        intervalMarkup += `:${strand}`;
      }
  
      return `<div>
        <div id="bed12-interval" style="display:block;">
          ${intervalMarkup}
        </div>
        <div id="bed12-score" style="display:block;">
          Score: ${realScore}
        </div>
        <div id="bed12-component" style="display:block;">
          ${itemRGBMarkup}
        </div>
        <div id="bed12-element-cartoon" style="display:block;">
          ${elementCartoon}
        </div>
      </div>`;
    }

    getMouseOverHtml(trackX, trackY) {
      if (!this.tilesetInfo) {
        return "";
      }

      const point = [trackX, trackY];

      for (const tile of this.visibleAndFetchedTiles()) {
        for (let i = 0; i < tile.allExonsForMouseOver.length; i++) {
          const rect = tile.allExonsForMouseOver[i][0];

          if (this.isPointInRectangle(rect, point)) {

            const transcript = tile.allExonsForMouseOver[i][1];
            
            switch (this.options.blockStyle) {
              case "directional":
              case "UCSC-like":
                {
                  let overlapType = "";
                  let overlapSubtype = null;
                  let overlapIndex = -1;
                  let overlapCount = -1;
                  const normPointX = this.getNormPointXWithinLocalRect(point[0], rect[0], rect[1]);

                  for (let testBlockIdx = 0; testBlockIdx < transcript.blocks.length; testBlockIdx++) {
                    const testBlock = transcript.blocks[testBlockIdx];
                    if ((normPointX >= testBlock.range[0]) && (normPointX <= testBlock.range[1])) {
                      overlapType = testBlock.type;
                      overlapIndex = testBlock.typeIdx;
                      switch (overlapType) {
                        case "Exon":
                          overlapCount = transcript.blockTypeCounts.exons;
                          overlapSubtype = (testBlock.subtype) ? testBlock.subtype : null;
                          break;
                        case "Intron":
                          overlapCount = transcript.blockTypeCounts.introns;
                          break;
                        case "5'UTR":
                        case "3'UTR":
                          overlapCount = transcript.blockTypeCounts.exons;
                          break;
                        default:
                          overlapCount = "";
                          break;
                      }
                      break;
                    }
                  }

                  function formattedOverlapByType(type, subtype, typeIdx, typeCount) {
                    let result = "";
                    switch (type) {
                      case "Exon":
                      case "Intron":
                      case "5'UTR":
                      case "3'UTR":
                        result = (subtype) ? `${type} (${subtype}): ${typeIdx} / ${typeCount}` : `${type}: ${typeIdx} / ${typeCount}`;
                        break
                      default:
                        break
                    }
                    return result;
                  }

                  const formattedOverlap = formattedOverlapByType(overlapType, overlapSubtype, overlapIndex, overlapCount);

                  if (this.areCodonsShown) {
                    for (let j = 0; j < tile.allCodonsForMouseOver.length; j++) {
                      const rectDim = tile.allCodonsForMouseOver[j][0];

                      if (this.isPointInRectangle(rectDim, point)) {

                        const aa = tile.allCodonsForMouseOver[j][1];

                        if (aa.key === "X") {
                          return `
                            <div style="text-align: center; padding:10px;">
                              <div style="text-transform: uppercase;">
                                <b>${aa.name}</b>
                              </div>
                              <div style="font-size: 10px;">
                                ${aa.codons.join(", ")}
                              </div>
                            </div>
                          `;
                        }

                        const essential = aa.essential
                          ? "essential"
                          : "non-essential";
                        return `
                          <div style="text-align: center; padding:10px;">
                            <div>
                              <img class="fit-picture" 
                                  src="${aa.image}" 
                                  style="width:125px;height:125px"
                                  alt="${aa.name}">
                            </div>
                            <div style="text-transform: uppercase; padding-top:6px;">
                              <b>${aa.name}</b>
                            </div>
                            <div style="padding-top:0px;">
                              <i>${aa.property}, ${essential}</i>
                            </div>
                            <div style="font-size: 10px;">
                              ${aa.codons.join(", ")}
                            </div>
                          </div>
                        `;
                      }
                    }
                  } 
                  else {
                    return `
                      <div>
                        <div><b>Transcript: ${transcript.transcriptName}</b></div>
                        <div>Position: ${transcript.chromName}:${transcript.txStart}-${transcript.txEnd}</div>
                        <div>Strand: ${transcript.strand}</div>
                        <div>${formattedOverlap}</div>
                      </div>
                    `;
                  }
                }
                break;
              case "boxplot":
                {
                  const bed12FieldsObj = {
                    chrom: transcript.chromName,
                    start: transcript.txStart,
                    end: transcript.txEnd,
                    id: transcript.rawTranscriptName,
                    score: transcript.importance,
                    strand: transcript.strand,
                    thickStart: transcript.startCodonPos,
                    thickEnd: transcript.stopCodonPos,
                    itemRGB: (transcript.itemRgb !== '.') ? transcript.itemRgb : '0,0,0',
                    blockCount: transcript.boxBlockCount,
                    blockSizes: transcript.boxBlockLengths,
                    blockStarts: transcript.boxBlockStarts,
                  };

                  const bed12Output = this.formattedBED12HTML(bed12FieldsObj);

                  if (bed12Output.length === 0) {
                    const hc = document.getElementsByClassName('higlass-main-content')[0];
                    if (hc) {
                      hc.style.cursor = 'grab';
                    }
                  }

                  return bed12Output;
                }
                break;
              default:
                throw new Error(
                  'Uncaught TypeError: Unknown blockStyle option (drawLabels, within tile data)'
                );
            }
          }
        }
      }

      return "";
    }

    exportSVG() {
      let track = null;
      let base = null;
  
      base = document.createElement('g');
      track = base;

      const clipPathId = slugid.nice();

      const gClipPath = document.createElement('g');
      gClipPath.setAttribute('style', `clip-path:url(#${clipPathId});`);

      track.appendChild(gClipPath);

      // define the clipping area as a polygon defined by the track's
      // dimensions on the canvas
      const clipPath = document.createElementNS(
        'http://www.w3.org/2000/svg',
        'clipPath'
      );
      clipPath.setAttribute('id', clipPathId);
      track.appendChild(clipPath);

      const clipPolygon = document.createElementNS(
        'http://www.w3.org/2000/svg',
        'polygon'
      );
      clipPath.appendChild(clipPolygon);

      clipPolygon.setAttribute(
        'points',
        `${this.position[0]},${this.position[1]} ` +
          `${this.position[0] + this.dimensions[0]},${this.position[1]} ` +
          `${this.position[0] + this.dimensions[0]},${this.position[1] +
            this.dimensions[1]} ` +
          `${this.position[0]},${this.position[1] + this.dimensions[1]} `
      );

      const output = document.createElement('g');
      
      output.setAttribute(
        'transform',
        `translate(${this.position[0]},${this.position[1]})`
      );
  
      gClipPath.appendChild(output);
      
      // paint a backing rectangle with options.backgroundColor
      const backingRect = document.createElement('rect');
      backingRect.setAttribute('fill', this.options.backgroundColor);
      backingRect.setAttribute('fill-opacity', '1');
      backingRect.setAttribute('stroke-opacity', '0');
      backingRect.setAttribute('width', this.dimensions[0]);
      backingRect.setAttribute('height', this.dimensions[1]);
      backingRect.setAttribute('x', '0');
      backingRect.setAttribute('y', '0');
      output.appendChild(backingRect);

      // We need to draw the lower order rectangles first (middle line)
      const paintOrders = [0,1];

      paintOrders.forEach(paintOrder => {

        this.visibleAndFetchedTiles()
          .filter(tile => tile.allExonsForSVG)
          .forEach(tile => {
            const gTile = document.createElement('g');
            gTile.setAttribute(
              'transform',
              `translate(${tile.rectGraphics.position.x},
              ${tile.rectGraphics.position.y})
              scale(${tile.rectGraphics.scale.x},
              ${tile.rectGraphics.scale.y})`
            );
    
            tile.allExonsForSVG
            .filter(rect => rect.paintOrder === paintOrder)
            .forEach(rect => {
              const r = document.createElement('path');
              const poly = rect.rect;

              let d = `M ${poly[0]} ${poly[1]}`;
              for (let i = 2; i < poly.length; i += 2) {
                d += ` L ${poly[i]} ${poly[i + 1]}`;
              }
              d += " Z";

              r.setAttribute('d', d);
              r.setAttribute('fill', rect.color);
              r.setAttribute('opacity', '1');
              gTile.appendChild(r);
            });
            output.appendChild(gTile);
          });

      });
      
      // We don't want to draw texts twice
      const alreadyDrawnTexts = [];

      const polyMargin = {
        'left' : 2,
        'right' : 2,
        'top' : 1,
        'bottom' : 1
      };
      const polyTriangleWidth = 5;
      const labelYOffsetAdjustment = 1;

      this.allTexts
        .filter(text => text.text.visible)
        .forEach(text => {

          if (alreadyDrawnTexts.includes(text.text.text)){
            return;
          }

          const poly = document.createElement('polygon');
          switch (this.options.blockStyle) {
            case "directional":
            case "UCSC-like": {
              poly.setAttribute('stroke-width', '1');
              poly.setAttribute('stroke-opacity', '1');
              // strand-orient the points
              const polyTopLeft = `${-polyMargin.left},${polyMargin.top}`;
              const polyTopRight = `${text.text.width + polyMargin.right},${polyMargin.top}`;
              const polyBottomRight = `${text.text.width + polyMargin.right},${-this.transcriptHeight*0.75 - polyMargin.bottom}`;
              const polyBottomLeft = `${-polyMargin.left},${-this.transcriptHeight*0.75 - polyMargin.bottom}`;
              let polyFillColor = (typeof this.options.labelBackgroundPlusStrandColor !== 'undefined') ? this.options.labelBackgroundPlusStrandColor : TranscriptsTrack.config.defaultOptions.labelBackgroundPlusStrandColor;
              switch (text.strand) {
                case '-':
                  const polyStrokeMinusColor = (typeof this.options.labelStrokeMinusStrandColor !== 'undefined') ? this.options.labelStrokeMinusStrandColor : TranscriptsTrack.config.defaultOptions.labelStrokeMinusStrandColor;
                  poly.setAttribute('stroke', polyStrokeMinusColor);
                  polyFillColor = (typeof this.options.labelBackgroundMinusStrandColor !== 'undefined') ? this.options.labelBackgroundMinusStrandColor : TranscriptsTrack.config.defaultOptions.labelBackgroundMinusStrandColor;
                  
                  const polyMiddleLeft = `${-polyMargin.left - polyTriangleWidth},${-this.transcriptHeight*0.75/2.0}`;
                  poly.setAttribute('points', `${polyMiddleLeft} ${polyTopLeft} ${polyTopRight} ${polyBottomRight} ${polyBottomLeft}`);
                  break;
                case '+':
                default:
                  const polyStrokePlusColor = (typeof this.options.labelStrokeMinusStrandColor !== 'undefined') ? this.options.labelStrokeMinusStrandColor : TranscriptsTrack.config.defaultOptions.labelStrokeMinusStrandColor;
                  poly.setAttribute('stroke', polyStrokePlusColor);
                  const polyMiddleRight = `${text.text.width + polyMargin.right + polyTriangleWidth},${-this.transcriptHeight*0.75/2.0}`;
                  poly.setAttribute('points', `${polyTopLeft} ${polyTopRight} ${polyMiddleRight} ${polyBottomRight} ${polyBottomLeft}`);
                  break;
              }
              switch (this.options.highlightTranscriptType) {
                case "none":
                  break;
                case "longestIsoform":
                  polyFillColor = (text.isLongestIsoform) ? this.options.highlightTranscriptLabelBackgroundColor : polyFillColor;
                  break;
                case "apprisPrincipalIsoform":
                  polyFillColor = (text.isApprisPrincipalIsoform) ? this.options.highlightTranscriptLabelBackgroundColor : polyFillColor;
                  break;
                default:
                  throw new Error(
                    'Uncaught TypeError: Unknown highlightTranscriptType option (transcript label background, SVG export)'
                  );
              }
              poly.setAttribute('fill', polyFillColor);
              break;
            }
            case "boxplot": {
              break;
            }
            default: {
              throw new Error(
                'Uncaught TypeError: Unknown blockStyle option (exportSVG, poly-attributes)'
              );
            }
          }

          const t = document.createElement('text');
          t.setAttribute('text-anchor', 'start');
          t.setAttribute('font-family', this.options.fontFamily);
          t.setAttribute('font-weight', text.text.style.fontWeight); // pull weight from entity
          t.setAttribute('font-size', `${this.options.labelFontSize}px`);
          t.setAttribute('fill', this.options.labelFontColor);
          t.innerHTML = text.text.text;
          alreadyDrawnTexts.push(text.text.text);

          const g = document.createElement('g');
          g.setAttribute('transform', `scale(${text.text.scale.x},1)`);
          switch (this.options.blockStyle) {
            case "directional":
            case "UCSC-like": {
              g.appendChild(poly);
              break;
            }
            case "boxplot": {
              break;
            }
            default: {
              throw new Error(
                'Uncaught TypeError: Unknown blockStyle option (exportSVG, poly-append)'
              );
            }
          }
          g.appendChild(t);
          g.setAttribute(
            'transform',
            `translate(${text.text.x},${text.text.y + this.transcriptHeight/4 + labelYOffsetAdjustment})scale(${text.text.scale.x},1)`
          );

          output.appendChild(g);
        });  
  
      return [base, base];
    }

  }
  return new TranscriptsTrackClass(...args);
};

const icon =
  '<svg width="20" height="20" xmlns="http://www.w3.org/2000/svg"> <g> <title>background</title> <rect fill="#fff" id="canvas_background" height="3.24996" width="3.24996" y="-1" x="-1"/> <g display="none" id="canvasGrid"> <rect fill="url(#gridpattern)" stroke-width="0" y="0" x="0" height="100%" width="100%" id="svg_2"/> </g> </g> <g> <title>Layer 1</title> <rect id="svg_1" height="20" width="20" stroke-width="0" stroke="#000" fill="#C0EAAF"/> <path stroke="#bdbfff" id="svg_4" d="m11.42795,10.0119l-4.10746,-6.99997l2.70509,0l4.10746,6.99997l-4.10746,6.99997l-2.70509,0l4.10746,-6.99997z" stroke-width="6" fill="#bdbfff"/> </g> </svg>';

// default
TranscriptsTrack.config = {
  type: "horizontal-transcripts",
  datatype: ["gene-annotation"],
  local: false,
  orientation: "1d-horizontal",
  thumbnail: new DOMParser().parseFromString(icon, "text/xml").documentElement,
  availableOptions: [
    "fontSize",
    "fontFamily",
    "transcriptSpacing",
    "transcriptHeight",
    "maxTexts",
    "maxRows",
    "plusStrandColor",
    "minusStrandColor",
    "utrColor",
    "labelBackgroundPlusStrandColor",
    "labelBackgroundMinusStrandColor",
    "labelFontColor",
    "labelFontSize",
    "labelFontWeight",
    "labelStrokePlusStrandColor",
    "labelStrokeMinusStrandColor",
    "startCollapsed",
    "showToggleTranscriptsButton",
    "trackHeightAdjustment",
    "sequenceData",
    "backgroundColor",
    "blockStyle",
    "highlightTranscriptType",
    "highlightTranscriptTrackBackgroundColor",
    "highlightTranscriptLabelBackgroundColor",
    "highlightTranscriptLabelFontWeight",
    "showHighlightedTranscriptsOnly",
    "isVisible",
    "allowResizeParentDiv",
  ],
  defaultOptions: {
    fontSize: 9,
    fontFamily: "Helvetica",
    transcriptSpacing: 2,
    transcriptHeight: 11,
    maxTexts: 100,
    maxRows: null,
    plusStrandColor: "#bdbfff",
    minusStrandColor: "#fabec2",
    utrColor: "#C0EAAF",
    labelBackgroundPlusStrandColor: "#ffffff",
    labelBackgroundMinusStrandColor: "#ffffff",
    labelFontColor: "#333333",
    labelFontSize: 10,
    labelFontWeight: "300",
    labelStrokePlusStrandColor: "#999999",
    labelStrokeMinusStrandColor: "#999999",
    startCollapsed: false,
    trackHeightAdjustment: "automatic",
    showToggleTranscriptsButton: true,
    backgroundColor: "#ffffff",
    blockStyle: "directional", // "directional" | "UCSC-like" | "boxplot"
    highlightTranscriptType: "none", // "none" | "longestIsoform" | "apprisPrincipalIsoform"
    highlightTranscriptTrackBackgroundColor: "#f0f0f0",
    highlightTranscriptLabelBackgroundColor: "#f0f0f0",
    highlightTranscriptLabelFontWeight: "700",
    showHighlightedTranscriptsOnly: false,
    isVisible: true,
    allowResizeParentDiv: true,
  },
};

export default TranscriptsTrack;
