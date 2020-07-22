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

const TranscritpsTrack = (HGC, ...args) => {
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
  const DARKGREY_HEX = colorToHex("#999999");

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
    const { flipText, fontSize, fontFamily, maxTexts } = options;
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
      const transcriptName = tsFormatted.transcriptName;
      const transcriptId = tsFormatted.transcriptId;
      const strand = tsFormatted.strand;

      td["transcriptId"] = transcriptId;

      // don't draw texts for the latter entries in the tile
      if (i >= maxTexts) return;

      const text = new HGC.libraries.PIXI.Text(transcriptName, {
        fontSize: `${fontSize}px`,
        fontFamily,
        fill: track.colors["labelFont"],
      });
      text.interactive = true;

      if (flipText) text.scale.x = -1;

      text.anchor.x = 0;
      text.anchor.y = 0.5;
      text.visible = false;

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

    // Compute the offsets of each exon, so that we can get codons accross exons
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

    // Add start and stop codon to the exon list and distingush between UTR and coding region later
    if (isProteinCoding) {
      exonOffsetStarts.push(startCodonPos, stopCodonPos);
      exonOffsetEnds.push(startCodonPos, stopCodonPos);

      exonOffsetStarts.sort();
      exonOffsetEnds.sort();
    }

    const xStartPos = track._xScale(txStart);
    const xEndPos = track._xScale(txEnd);

    const width = xEndPos - xStartPos;
    const yMiddle = centerY;

    const polys = [];
    const polysSVG = []; // holds all polygons that need to be drawn for SVG export

    graphics.beginFill(track.colors.intron);
    // draw the middle line
    let poly = [
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

    // For mouseOver
    polys.push([xStartPos, xStartPos + width, topY, topY + height]);

    // For SVG export
    polysSVG.push({
      rect: poly,
      color: track.colors.intronHEX,
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

      const colorUsed = isNonCodingOrUtr ? track.colors.utr : track.colors[strand];
      const colorUsedSVG = isNonCodingOrUtr ? track.options.utrColor : track.colors[strand+"HEX"];
      
      graphics.beginFill(colorUsed);
      const xStart = track._xScale(exonStart);
      const localWidth = Math.max(
        2,
        track._xScale(exonEnd) - track._xScale(exonStart)
      );

      let minX = xStartPos;
      let maxX = xEndPos;
      let localPoly = null;
      let localRect = null; // without direction for mouseOver

      if (strand === "+") {
        const rectStartX = Math.min(xStart, maxX);
        const rectStartX2 = Math.max(rectStartX - 5, xStartPos);
        const rectEndX = Math.min(xStart + localWidth, maxX);
        const rectEndX2 = Math.max(rectEndX - 5, xStartPos);

        localPoly = [
          rectStartX,
          topY,
          rectEndX2,
          topY,
          rectEndX,
          topY + height / 2,
          rectEndX2,
          topY + height,
          rectStartX2,
          topY + height,
          rectStartX,
          topY + height / 2,
          rectStartX2,
          topY,
        ];

        localRect = [rectStartX, rectEndX, topY, topY + height];
      } else {
        const rectStartX = Math.max(xStart, minX);
        const rectStartX2 = Math.min(rectStartX + 5, xEndPos);
        const rectEndX = Math.max(xStart + localWidth, minX);
        const rectEndX2 = Math.min(rectEndX + 5, xEndPos);

        localPoly = [
          rectStartX,
          topY + height / 2,
          rectStartX2,
          topY,
          rectEndX2,
          topY,
          rectEndX,
          topY + height / 2,
          rectEndX2,
          topY + height,
          rectStartX2,
          topY + height,
          rectStartX,
          topY + height / 2,
        ];

        localRect = [rectStartX, rectEndX, topY, topY + height];
      }

      graphics.drawPolygon(localPoly);
      polys.push(localRect);

      // For SVG export
      polysSVG.push({
        rect: localPoly,
        color: colorUsedSVG,
        paintOrder: 1
      });
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


      if (track.areTranscriptsHidden && track.transcriptInfo[transcriptId].displayOrder !== 0){
        return;
      };

      let centerYOffset =
        track.transcriptInfo[transcriptId].displayOrder *
        (height + strandSpacing);

      if (track.options.showToggleTranscriptsButton) {
        centerYOffset += track.toggleButtonHeight;
      }

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
    });
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
    const width = track._xScale(tileX + tileWidth) - track._xScale(tileX);
    const height = track.dimensions[1];
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

  class TranscritpsTrackClass extends HGC.tracks
    .HorizontalGeneAnnotationsTrack {
    constructor(context, options) {
      super(context, options);
      const { animate } = context;

      this.trackId = this.id;

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

      if (this.options.sequenceData !== undefined) {
        this.sequenceLoader = new SequenceLoader(
          this.options.sequenceData.fastaUrl,
          this.options.sequenceData.faiUrl
        );
        if (!this.pixiTexts) {
          this.pixiTexts = initializePixiTexts(this.codonTextOptions, HGC);
        }
      }

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
      this.colors["labelBackground"] = colorToHex(
        this.options.labelBackgroundColor
      );
    }

    initTile(tile) {
      externalInitTile(this, tile, {
        flipText: this.flipText,
        fontSize: this.fontSize,
        fontFamily: this.options.fontFamily,
        maxTexts: this.options.maxTexts,
      });

      // We have to rerender everything since the vertical position
      // of the tracks might have changed accross tiles
      this.rerender(this.options, true);
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
      let height = 0;
      const trackMargin = 10;

      if (this.areTranscriptsHidden) {
        height = this.toggleButtonHeight + 
          Math.min(1, this.numTranscriptRows) * (this.transcriptHeight + this.transcriptSpacing) + 
          trackMargin;
      } else {
        const tbh = this.options.showToggleTranscriptsButton
          ? this.toggleButtonHeight
          : 0;

        height =
          this.numTranscriptRows *
            (this.transcriptHeight + this.transcriptSpacing) +
          tbh +
          trackMargin;
      }

      this.trackHeightOld = this.trackHeight;
      this.trackHeight = height;

      return height;
    }

    adjustTrackHeight() {
      this.computeTrackHeight();

      if (this.trackHeightOld === this.trackHeight) {
        return false;
      }

      this.pubSub.publish("trackDimensionsModified", {
        height: this.trackHeight,
        resizeParentDiv: true,
        trackId: this.trackId,
      });

      return true;
    }

    formatTranscriptData(ts) {
      const strand = ts[5];
      const stopCodonPos = ts[12] === "." ? "." : (strand === "+" ? +ts[12] + 2 : +ts[12] - 1);
      const startCodonPos = ts[11] === "." ? "." : (strand === "+" ? +ts[11] - 1 : +ts[11] + 2);
      const exonStarts = ts[9].split(",").map((x) => +x - 1);
      const exonEnds = ts[10].split(",").map((x) => +x);
      const txStart = +ts[1] - 1;
      const txEnd = +ts[2];

      const result = {
        transcriptId: this.transcriptId(ts),
        transcriptName: ts[3],
        txStart: txStart,
        txEnd: txEnd,
        strand: strand,
        chromName: ts[0],
        codingType: ts[8],
        exonStarts: exonStarts,
        exonEnds: exonEnds,
        startCodonPos: startCodonPos,
        stopCodonPos: stopCodonPos,
        importance: +ts[4],
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
        visibleTranscripts.push(visibleTranscriptsObj[tsId]);
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

          if (this.transcriptPositionInfo[dpo] === undefined) {
            this.transcriptPositionInfo[dpo] = [];
          }
          this.transcriptPositionInfo[dpo].push([+ts[1], +ts[2], ts[3]]);

          const tsFormatted = this.formatTranscriptData(ts);

          const tInfo = {
            transcriptId: tsFormatted.transcriptId,
            transcriptName: tsFormatted.transcriptName,
            txStart: tsFormatted.txStart,
            txEnd: tsFormatted.txEnd,
            strand: tsFormatted.strand,
            chromName: tsFormatted.chromName,
            codingType: tsFormatted.codingType,
            exonStarts: tsFormatted.exonStarts,
            exonEnds: tsFormatted.exonEnds,
            startCodonPos: tsFormatted.startCodonPos,
            stopCodonPos: tsFormatted.stopCodonPos,
            displayOrder: dpo,
            importance: tsFormatted.importance,
          };
          this.transcriptInfo[tInfo.transcriptId] = tInfo;
        });

      this.numTranscriptRows = Object.keys(this.transcriptPositionInfo).length;

      // Update the button text.
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


      const strandCenterY =
        this.transcriptHeight / 2 + this.transcriptSpacing / 2;

      const renderContext = [
        this,
        tile,
        strandCenterY,
        this.transcriptHeight,
        this.transcriptSpacing,
      ];

      renderTranscriptExons(plusTranscripts, ...renderContext);
      renderTranscriptExons(minusTranscripts, ...renderContext);

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

            if (this.areTranscriptsHidden && transcript.displayOrder !== 0) {
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

              const startPosX = this._xScale(codon.posStart + chrOffset);
              const endPosX = this._xScale(codon.posEnd + chrOffset);

              // if the codon is not visible, don't bother
              if (endPosX < 0 || startPosX > this.dimensions[0]) {
                continue;
              }

              //const availableSpace = this._xScale((codon.posStart)) - this._xScale((codon.posEnd + 1));
              //console.log(availableSpace);
              if (codonWidth < this.minCodonDistance + 10) {
                const xMiddle =
                  this._xScale(
                    (codon.posStart + codon.posEnd + 1) / 2 + chrOffset
                  ) -
                  codon.widthAbbrev / 2; //(codon.posStart + codon.posEnd + 1) / 2 + chrOffset
                codon.spriteAbbrev.position.x = xMiddle;
                codon.spriteAbbrev.position.y = yMiddle;
                graphics.addChild(codon.spriteAbbrev);
              } else {
                const xMiddle =
                  this._xScale(
                    (codon.posStart + codon.posEnd + 1) / 2 + chrOffset
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

            if (this.areTranscriptsHidden && transcript.displayOrder !== 0) {
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

              const startPosX = this._xScale(codon.posStart + chrOffset);
              const endPosX = this._xScale(codon.posEnd + chrOffset + 1);

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
          tile.labelBgGraphics.beginFill(
            typeof this.options.labelBackgroundColor !== "undefined"
              ? this.colors["labelBackground"]
              : WHITE_HEX
          );

          if (!tile.initialized) return;

          tile.tileData.forEach((td) => {
            // tile probably hasn't been initialized yet
            if (!tile.texts) return;

            const transcriptId = td.transcriptId;

            if (this.transcriptInfo[transcriptId] === undefined) return;

            const transcript = this.transcriptInfo[transcriptId];
            const text = tile.texts[transcriptId];

            if (!text) return;

            if (this.areTranscriptsHidden && transcript.displayOrder !== 0) {
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

            // take care of label positioning at start or end of transcripts
            if (transcript.strand === "+") {
              text.position.x = Math.max(
                this._xScale(this.xScale().domain()[0]) + TEXT_MARGIN,
                this._xScale(txStart) -
                  tile.textWidths[transcriptId] -
                  2 * TEXT_MARGIN -
                  CARET_MARGIN
              );
            } else {
              text.position.x = Math.max(
                this._xScale(this.xScale().domain()[0]) +
                  TEXT_MARGIN +
                  CARET_WIDTH,
                this._xScale(txStart) -
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
              this._xScale(txEnd) - marginRight
            );

            text.position.y = textYMiddle;

            // Determine if the current text should be hidden
            let showText = true;
            const dpo = transcript.displayOrder;

            this.transcriptPositionInfo[dpo]
              .filter((ts) => {
                // Check the ones that are left of the current transcript
                return ts[1] < transcript.txStart;
              })
              .forEach((ts) => {
                const endOfTranscript = this._xScale(ts[1] + chrOffset);

                if (endOfTranscript > text.position.x - 4 * TEXT_MARGIN) {
                  showText = false;
                }
              });

            if (showText) {
              text.visible = true;

              this.allBoxes.push([
                text.position.x - TEXT_MARGIN,
                textYMiddle,
                tile.textWidths[transcriptId] + 2 * TEXT_MARGIN,
                this.transcriptHeight,
                transcript.transcriptName,
              ]);

              this.allTexts.push({
                importance: transcript.importance,
                text,
                caption: transcript.transcriptName,
                strand: transcript.strand,
              });

              allTiles.push(tile.labelBgGraphics);
            } else {
              text.visible = false;
            }
          });
        });

      this.renderTextBg(this.allBoxes, this.allTexts, allTiles);
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
            allTiles[i].beginFill(DARKGREY_HEX);

            // Directional label
            let polyBorder = [xl, yb, xr, yb, xr + 5, minY, xr, yt, xl, yt];

            allTiles[i].drawPolygon(polyBorder);

            allTiles[i].beginFill(this.colors["labelBackground"]);

            // Directional label
            let poly = [xlm, ybm, xr, ybm, xrm + 3, minY, xr, ytm, xlm, ytm];
            allTiles[i].drawPolygon(poly);
          } else {
            allTiles[i].beginFill(WHITE_HEX);

            let polyBg = [xl - 5, yb, xr, yb, xr, yt, xl - 5, yt];

            allTiles[i].drawPolygon(polyBg);

            allTiles[i].beginFill(DARKGREY_HEX);
            // Directional label
            let polyBorder = [xl - 5, minY, xl, yb, xr, yb, xr, yt, minX, yt];

            allTiles[i].drawPolygon(polyBorder);

            allTiles[i].beginFill(this.colors["labelBackground"]);

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
            } else {
              return `
                <div>
                  <div><b>Transcript: ${transcript.transcriptName}</b></div>
                  <div>Position: ${transcript.chromName}:${transcript.txStart}-${transcript.txEnd}</div>
                  <div>Strand: ${transcript.strand}</div>
                </div>
              `;
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

              r.setAttribute('d', d);
              r.setAttribute('fill', rect.color);
              r.setAttribute('opacity', '1');
              gTile.appendChild(r);
            });
            output.appendChild(gTile);
          });

      });
      
      // We dont want to draw textsw twice
      const allreadyDrawnTexts = [];

      this.allTexts
        .filter(text => text.text.visible)
        .forEach(text => {

          if(allreadyDrawnTexts.includes(text.text.text)){
            return;
          }

          const g = document.createElement('g');
          const t = document.createElement('text');
          t.setAttribute('text-anchor', 'start');
          t.setAttribute('font-family', this.options.fontFamily);
          t.setAttribute('font-size', `${this.fontSize}px`);
  
          g.setAttribute('transform', `scale(${text.text.scale.x},1)`);
  
          t.setAttribute('fill', this.colors["labelFont"]);
          t.innerHTML = text.text.text;
          allreadyDrawnTexts.push(text.text.text);
  
          g.appendChild(t);
          g.setAttribute(
            'transform',
            `translate(${text.text.x},${text.text.y + this.transcriptHeight/4})scale(${text.text.scale.x},1)`
          );
          output.appendChild(g);
        });
  
      return [base, base];
    }

  }
  return new TranscritpsTrackClass(...args);
};

const icon =
  '<svg width="20" height="20" xmlns="http://www.w3.org/2000/svg"> <g> <title>background</title> <rect fill="#fff" id="canvas_background" height="3.24996" width="3.24996" y="-1" x="-1"/> <g display="none" id="canvasGrid"> <rect fill="url(#gridpattern)" stroke-width="0" y="0" x="0" height="100%" width="100%" id="svg_2"/> </g> </g> <g> <title>Layer 1</title> <rect id="svg_1" height="20" width="20" stroke-width="0" stroke="#000" fill="#C0EAAF"/> <path stroke="#bdbfff" id="svg_4" d="m11.42795,10.0119l-4.10746,-6.99997l2.70509,0l4.10746,6.99997l-4.10746,6.99997l-2.70509,0l4.10746,-6.99997z" stroke-width="6" fill="#bdbfff"/> </g> </svg>';

// default
TranscritpsTrack.config = {
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
    "plusStrandColor",
    "minusStrandColor",
    "utrColor",
    "labelBackgroundColor",
    "labelFontColor",
    "startCollapsed",
    "showToggleTranscriptsButton",
    "trackHeightAdjustment",
    "sequenceData",
  ],
  defaultOptions: {
    fontSize: 9,
    fontFamily: "Helvetica",
    transcriptSpacing: 2,
    transcriptHeight: 11,
    maxTexts: 20,
    plusStrandColor: "#bdbfff",
    minusStrandColor: "#fabec2",
    utrColor: "#C0EAAF",
    labelBackgroundColor: "#ffffff",
    labelFontColor: "#333333",
    startCollapsed: true,
    trackHeightAdjustment: "automatic",
    showToggleTranscriptsButton: true,
  },
};

export default TranscritpsTrack;
