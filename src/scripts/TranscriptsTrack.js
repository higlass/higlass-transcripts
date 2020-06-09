//import { scaleLinear, scaleOrdinal, schemeCategory10 } from "d3-scale";
//import { color } from "d3-color";
import boxIntersect from "box-intersect";
import classifyPoint from "robust-point-in-polygon";
import { AMINO_ACIDS, CODONS } from './configs';
import { TextStyle } from "pixi.js";
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
  const { colorToHex, trackUtils } = HGC.utils;

  // these are default values that are overwritten by the track's options

  const MAX_TEXTS = 100;
  const WHITE_HEX = colorToHex("#ffffff");

  const MAX_GENE_ENTRIES = 50;
  const MAX_FILLER_ENTRIES = 5000;

  /**
   * Initialize a tile. Pulled out from the track so that it
   * can be modified without having to modify the track
   * object (e.g. in an Observable notebooke)
   *
   * @param  {HorizontalGeneAnnotationsTrack} track   The track object
   * @param  {Object} tile    The tile to render
   * @param  {Object} options The track's options
   */
  function externalInitTile(track, tile, options) {
    const {
      flipText,
      fontSize,
      fontFamily,
      maxGeneEntries,
      maxFillerEntries,
      maxTexts,
    } = options;
    // create texts
    tile.texts = {};

    tile.rectGraphics = new HGC.libraries.PIXI.Graphics();
    tile.rectMaskGraphics = new HGC.libraries.PIXI.Graphics();
    tile.textBgGraphics = new HGC.libraries.PIXI.Graphics();
    tile.textGraphics = new HGC.libraries.PIXI.Graphics();

    tile.graphics.addChild(tile.rectGraphics);
    tile.graphics.addChild(tile.rectMaskGraphics);
    tile.graphics.addChild(tile.textBgGraphics);
    tile.graphics.addChild(tile.textGraphics);

    tile.rectGraphics.mask = tile.rectMaskGraphics;

    if (!tile.tileData.sort) return;

    tile.tileData.sort((a, b) => b.importance - a.importance);

    // const geneEntries = tile.tileData
    //   .filter(td => td.type !== 'filler')
    //   .slice(0, maxGeneEntries);

    //tile.tileData = geneEntries;

    tile.tileData.forEach((td, i) => {
      const geneInfo = td.fields;
      const geneName = geneInfo[3];
      const geneId = track.transcriptId(geneInfo);
      const strand = td.strand || geneInfo[5];

      tile.textWidths = {};
      tile.textHeights = {};

      // don't draw texts for the latter entries in the tile
      if (i >= maxTexts) return;

      const text = new HGC.libraries.PIXI.Text(geneName, {
        fontSize: `${fontSize}px`,
        fontFamily,
        fill: track.colors["labelFont"],
      });
      text.interactive = true;

      if (flipText) text.scale.x = -1;

      text.anchor.x = 0;
      text.anchor.y = 0.5;
      text.visible = false;

      tile.texts[geneId] = text; // index by geneName
      tile.texts[geneId].strand = strand;
      tile.textGraphics.addChild(text);
    });

    tile.initialized = true;
  }

  /** Draw the exons within a gene */
  function drawExons(
    track,
    graphics,
    txStart,
    txEnd,
    geneInfo,
    chrOffset,
    centerY,
    height,
    strand,
    alpha
  ) {
    const topY = centerY - height / 2;

    const exonStarts = geneInfo[12];
    const exonEnds = geneInfo[13];
    const isProteinCoding = geneInfo[8] === "protein_coding";
    const startCodonPos = isProteinCoding ? +geneInfo[14] : -1;
    const stopCodonPos = isProteinCoding ? +geneInfo[15] : -1;

    let exonOffsetStarts = exonStarts.split(",").map((x) => +x + chrOffset);
    let exonOffsetEnds = exonEnds.split(",").map((x) => +x + chrOffset);

    // Add start and stop codon to the exon list and distingush between UTR and coding reagion later
    if(isProteinCoding){
      exonOffsetStarts.push(startCodonPos, stopCodonPos);
      exonOffsetEnds.push(startCodonPos, stopCodonPos);
      //exonOffsetStarts = [... new Set(exonOffsetStarts.sort())];
      exonOffsetStarts.sort()
      //exonOffsetEnds = [... new Set(exonOffsetEnds.sort())];
      exonOffsetEnds.sort()
    }
    //console.log(exonOffsetStarts);
    //console.log(exonOffsetEnds);

    const xStartPos = track._xScale(txStart);
    const xEndPos = track._xScale(txEnd);

    const width = xEndPos - xStartPos;
    const yMiddle = centerY;

    const polys = [];

    graphics.beginFill(track.colors.intron, alpha);
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
    polys.push(poly);

    // the distance between the mini-triangles
    const triangleInterval = 3 * height;

    // the first triangle (arrowhead) will be drawn in renderGeneSymbols
    for (
      let j = Math.max(track.position[0], xStartPos) + triangleInterval;
      j < Math.min(track.position[0] + track.dimensions[0], xStartPos + width);
      j += triangleInterval
    ) {
      if (strand === "+") {
        poly = [
          j,
          yMiddle - track.miniTriangleHeight,
          j + track.miniTriangleHeight,
          yMiddle,
          j,
          yMiddle + track.miniTriangleHeight,
        ];
      } else {
        poly = [
          j,
          yMiddle - track.miniTriangleHeight,
          j - track.miniTriangleHeight,
          yMiddle,
          j,
          yMiddle + track.miniTriangleHeight,
        ];
      }

      polys.push(poly);
      graphics.drawPolygon(poly);
    }

    // draw the actual exons
    for (let j = 0; j < exonOffsetStarts.length; j++) {
      const exonStart = exonOffsetStarts[j];
      const exonEnd = exonOffsetEnds[j];

      const isNonCodingOrUtr = 
        !isProteinCoding ||
        (strand === "+" && (
        exonEnd <= startCodonPos ||
        exonStart >= stopCodonPos
        )) ||
        (strand === "-" && (
          exonStart >= startCodonPos ||
          exonEnd <= stopCodonPos
          ))
      //console.log(isNonCodingOrUtr);

      if(isNonCodingOrUtr){
        graphics.beginFill(track.colors.utr, alpha);
      }
      else{
        graphics.beginFill(track.colors[strand], alpha);
      }
      


      const xStart = track._xScale(exonStart);
      const localWidth = Math.max(
        1,
        track._xScale(exonEnd) - track._xScale(exonStart)
      );

      // we're not going to draw rectangles over the arrowhead
      // at the start of the gene
      let minX = xStartPos;
      let maxX = xEndPos;
      let localPoly = null;

      if (strand === "+") {
        maxX = xEndPos - track.transcriptHHeight;
        localPoly = [
          Math.min(xStart, maxX),
          topY,
          Math.min(xStart + localWidth, maxX),
          topY,
          Math.min(xStart + localWidth, maxX),
          topY + height,
          Math.min(xStart, maxX),
          topY + height,
          Math.min(xStart, maxX),
          topY,
        ];
      } else {
        minX = xStartPos + track.transcriptHHeight;
        localPoly = [
          Math.max(xStart, minX),
          topY,
          Math.max(xStart + localWidth, minX),
          topY,
          Math.max(xStart + localWidth, minX),
          topY + height,
          Math.max(xStart, minX),
          topY + height,
          Math.max(xStart, minX),
          topY,
        ];
      }

      polys.push(localPoly);
      graphics.drawPolygon(localPoly);
     
    }

    return polys;
  }

  

  /** Draw the arrowheads at the ends of genes */
  function renderGeneSymbols(
    genes,
    track,
    tile,
    oldGraphics,
    xScale,
    color,
    alpha,
    centerY,
    height,
    strandSpacing
  ) {
    genes.forEach((gene) => {
      const transcriptId = track.transcriptId(gene.fields);
      let centerYOffset =
        track.transcriptInfo[transcriptId].displayOrder *
        (height + strandSpacing);

      if (track.options.showToggleTranscriptsButton) {
        centerYOffset += track.toggleButtonHeight;
      }

      const topY = centerY + centerYOffset - height / 2;
      const xStart = track._xScale(gene.xStart);
      const xEnd = track._xScale(gene.xEnd);

      const graphics = new HGC.libraries.PIXI.Graphics();
      tile.rectGraphics.addChild(graphics);

      graphics.beginFill(color, alpha);
      //graphics.beginFill(color, 1.0);
      graphics.interactive = true;
      graphics.buttonMode = true;
      //graphics.mouseup = (evt) => geneClickFunc(evt, track, gene);

      let poly = [];
      if (gene.strand === "+" || gene.fields[5] === "+") {
        const pointerStart = Math.max(xStart, xEnd - track.transcriptHHeight);
        const pointerEnd = pointerStart + track.transcriptHHeight;

        poly = [
          pointerStart,
          topY,
          pointerEnd,
          topY + track.transcriptHHeight,
          pointerStart,
          topY + track.transcriptHeight,
        ];
      } else {
        const pointerStart = Math.min(xEnd, xStart + track.transcriptHHeight);
        const pointerEnd = pointerStart - track.transcriptHHeight;

        poly = [
          pointerStart,
          topY,
          pointerEnd,
          topY + track.transcriptHHeight,
          pointerStart,
          topY + track.transcriptHeight,
        ];
      }

      graphics.drawPolygon(poly);
      tile.allRects.push([poly, gene.strand, gene]);
      //centerYOffset = centerYOffset + height + strandSpacing;
    });
  }

  function renderGeneExons(
    genes,
    track,
    tile,
    rectGraphics,
    xScale,
    color,
    alpha,
    centerY,
    height,
    strandSpacing
  ) {
    genes.forEach((gene) => {
      const geneInfo = gene.fields;
      const chrOffset = +gene.chrOffset;

      const transcriptId = track.transcriptId(geneInfo);
      let centerYOffset =
        track.transcriptInfo[transcriptId].displayOrder *
        (height + strandSpacing);

      if (track.options.showToggleTranscriptsButton) {
        centerYOffset += track.toggleButtonHeight;
      }

      
      const graphics = new HGC.libraries.PIXI.Graphics();
      tile.rectGraphics.addChild(graphics);

      graphics.beginFill(color, alpha);
      graphics.interactive = true;
      graphics.buttonMode = true;
      //graphics.mouseup = (evt) => geneClickFunc(evt, track, gene);

      tile.allRects = tile.allRects.concat(
        drawExons(
          track,
          graphics,
          gene.xStart,
          gene.xEnd,
          geneInfo,
          chrOffset, // not used for now because we have just one chromosome
          centerY + centerYOffset,
          height,
          gene.strand || gene.fields[5],
          alpha
        ).map((x) => [x, gene.strand, gene])
      );

      //centerYOffset = centerYOffset + height + strandSpacing;
    });
  }

  function renderGenes(
    genes,
    track,
    tile,
    graphics,
    xScale,
    color,
    alpha,
    centerY,
    height,
    strandSpacing
  ) {
    renderGeneSymbols(
      genes,
      track,
      tile,
      graphics,
      xScale,
      color,
      alpha,
      centerY,
      height,
      strandSpacing
    );
    renderGeneExons(
      genes,
      track,
      tile,
      graphics,
      xScale,
      color,
      alpha,
      centerY,
      height,
      strandSpacing
    );
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
      track.pToggleButton.children[0].alpha = 0.8;
      document.body.style.cursor = "pointer"; // I guess that's not very elegant
    } else if (overOrOut === "out") {
      track.pToggleButton.children[0].alpha = 0.5;
      document.body.style.cursor = "default";
    }
    requestAnimationFrame(track.animate);
  };

  const toggleBtnClick = (event, track) => {
    //console.log("click");
    //const text = track.pToggleButton.children[0];
    //console.log(text)
    if (!track.areTranscriptsHidden) {
      track.pToggleButton.children[0].text = "SHOW TRANSCRIPTS";
      track.areTranscriptsHidden = true;
    } else {
      track.pToggleButton.children[0].text = "HIDE TRANSCRIPTS";
      track.areTranscriptsHidden = false;
    }

    track.pubSub.publish("trackDimensionsModified", {
      height: track.computeTrackHeight(),
      resizeParentDiv: true,
      trackId: track.trackId,
    });

  };

  /** Create a preventing this track from drawing outside of its
   * visible area
   */
  function renderToggleBtn(track) {
    if (
      !track.options.showToggleTranscriptsButton ||
      track.hasToggleBtnBeenRendered
    ) {
      return;
    }

    track.pToggleButton = new HGC.libraries.PIXI.Graphics();
    track.pToggleButton.interactive = true;
    track.pToggleButton.buttonMode = true;
    track.pForeground.removeChildren();
    track.pForeground.addChild(track.pToggleButton);
    track.pToggleButton.clear();
    track.pToggleButton.removeChildren();
    //console.log(track.pForeground);

    //track.pToggleButton.beginFill(colorToHex('#ffffff'), 1);

    const text = new HGC.libraries.PIXI.Text("HIDE TRANSCRIPTS", {
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

    track.pToggleButton.addChild(text);
    track.hasToggleBtnBeenRendered = true;
    //console.log(track.pToggleButton);
    //track.pToggleButton.drawRect(0, 0, track.dimensions[0], track.toggleButtonHeight);
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

      this.areTranscriptsHidden = false;

      if(this.options.sequenceData !== undefined){
        this.sequenceLoader = new SequenceLoader(
          this.options.sequenceData.fastaUrl,
          this.options.sequenceData.faiUrl)
      }

      this.transcriptInfo = {};
      this.exonInformation = {};

      console.log(this.sequenceLoader);
      console.log("construct");
      
      //console.log(context);
      //console.log(this);
    }

    initOptions(){

      this.fontSize = +this.options.fontSize;
      this.transcriptHeight = +this.options.transcriptHeight ;

      this.transcriptSpacing = +this.options.transcriptSpacing;
      this.geneStrandHSpacing = this.transcriptSpacing / 2;
      this.transcriptHHeight = this.transcriptHeight / 2;

      this.miniTriangleHeight = 4 * (this.transcriptHeight) / 9;

      this.toggleButtonHeight = 26;

      this.colors = {};
      this.colors["+"] = colorToHex(this.options.plusStrandColor);
      this.colors["-"] = colorToHex(this.options.minusStrandColor);
      this.colors["utr"] = colorToHex(this.options.utrColor);
      this.colors["labelFont"] = colorToHex(this.options.labelFontColor);
      this.colors["black"] = colorToHex("#000000");
      this.colors["intron"] = colorToHex("#CFCFCF");
      this.colors["labelBackground"] = colorToHex(this.options.labelBackgroundColor);


    }

    initTile(tile) {
      externalInitTile(this, tile, {
        flipText: this.flipText,
        fontSize: this.fontSize,
        fontFamily: this.options.fontFamily,
        maxGeneEntries: MAX_GENE_ENTRIES,
        maxFillerEntries: MAX_FILLER_ENTRIES,
        maxTexts: MAX_TEXTS,
      });

      console.log("init");
      console.log(tile);
      //console.log(this.tilesetInfo);
      //console.log(this.zoomLevel);
      //console.log(tile.tileId.split(".")[1]);
      

      // We have to rerender everything since the vertical position
      // of the tracks might have changed accross tiles
      this.rerender(this.options, true);
      //this.renderTile(tile);
    }

    updateExonInformation(){

      if(this.zoomLevel !== this.tilesetInfo.max_zoom || this.sequenceLoader === undefined){
        return
      }

      this.exonInformation = {};

      this.visibleAndFetchedTiles()
        // tile hasn't been drawn properly because we likely got some
        // bogus data from the server
        .filter(tile => tile.drawnAtScale)
        .forEach((tile) => {
          tile.tileData.forEach((ts) => {
            const tsId = this.transcriptId(ts.fields);
            visibleTranscriptsObj[tsId] = ts.fields;
          });
        });


      if(this.zoomLevel === this.tilesetInfo.max_zoom && this.sequenceLoader !== undefined){
        const tileId = +tile.tileId.split(".")[1];
        this.sequenceLoader
          .getTile(this.zoomLevel, tileId, this.tilesetInfo)
          .then((values) => {

            
            
            const sequence =  values[0];


            //console.log(sequence);

            const tileWidth = +this.tilesetInfo.max_width / 2 ** +this.zoomLevel ;

            // get the bounds of the tile
            const tileStart = this.tilesetInfo.min_pos[0] + tileId * tileWidth;// computed too many times - improve
            const tileEnd = this.tilesetInfo.min_pos[0] + (tileId+1) * tileWidth;// computed too many times - improve
            console.log(tileStart);

            const visibleExons = []
            tile.tileData.forEach(td => {
              const exonStarts = td.fields[12].split(",").map((x) => +x);
              const exonEnds = td.fields[13].split(",").map((x) => +x);
              
              const entry = {
                transcriptName: td.fields[3],
                transcriptId: td.fields[7],
                exonStart: [],
                exonEnd: [],
              };
              for(let i = 0; i < exonStarts.length; i++){
                if(
                  (exonStarts[i] <= tileStart && exonEnds[i] >= tileStart) ||
                  (exonStarts[i] >= tileStart && exonEnds[i] <= tileEnd) ||
                  (exonStarts[i] <= tileEnd && exonEnds[i] >= tileEnd)
                  ){
                    entry.exonStart.push(exonStarts[i]);
                    entry.exonEnd.push(exonEnds[i]);
                }
              }
              visibleExons.push(entry);
              
            });
            
            console.log("Visible exonds", visibleExons);

          });
        
      }
        
    }

    /** cleanup */
    destroyTile(tile) {
      tile.rectGraphics.destroy();
      tile.rectMaskGraphics.destroy();
      tile.textGraphics.destroy();
      tile.textBgGraphics.destroy();
      tile.graphics.destroy();
    }

    computeTrackHeight() {
      let height = 0;
      const trackMargin = 10;

      if(this.areTranscriptsHidden){
        height = this.toggleButtonHeight + trackMargin
      }
      else{
        const tbh = this.options.showToggleTranscriptsButton
        ? this.toggleButtonHeight
        : 0;

        height =
        this.numTranscriptRows * (this.transcriptHeight + this.transcriptSpacing) +
        tbh +
        trackMargin;
      }

      this.trackHeightOld = this.trackHeight;
      this.trackHeight = height;

      return height;
    }

    adjustTrackHeight() {

      this.computeTrackHeight()
      //console.log(this.trackHeightOld, this.trackHeight);
      if(this.trackHeightOld === this.trackHeight){
        return false
      };

      this.pubSub.publish("trackDimensionsModified", {
        height: this.trackHeight,
        resizeParentDiv: true,
        trackId: this.trackId,
      });

      return true
  
      //requestAnimationFrame(track.animate);
    };

    updateTranscriptInfo() {
      // get all visible transcripts
      const visibleTranscriptsObj = {};

      this.visibleAndFetchedTiles()
        // tile hasn't been drawn properly because we likely got some
        // bogus data from the server
        //.filter(tile => tile.drawnAtScale)
        .forEach((tile) => {
          tile.tileData.forEach((ts) => {
            const tsId = this.transcriptId(ts.fields);
            visibleTranscriptsObj[tsId] = ts.fields;
          });
        });

      const visibleTranscripts = [];
      for (const tsId in visibleTranscriptsObj) {
        visibleTranscripts.push(visibleTranscriptsObj[tsId]);
      }

      // Delete transcripts that are not visible anymore and collext display orders
      // for(const tsId in this.transcriptInfo) {
      //   if(visibleTranscripts[tsId] === undefined){
      //     delete this.transcriptInfo[tsId];
      //   }
      // }
      this.transcriptInfo = {};
      this.transcriptPositionInfo = {};

      // const numTranscripts = Object.keys(visibleTranscripts).length;
      // const availableDisplayOrders = [...Array(numTranscripts).keys()]
      //   .filter(x => !displayOrders.includes(x));

      let displayOrder = 0;
      this.numTranscriptRows = 0;
      visibleTranscripts
        .sort(function (a, b) {
          return +a[1] - b[1];
        })
        .forEach((ts) => {

          const dpo = this.calculateTranscriptRowNumber(this.transcriptPositionInfo,+ts[1],+ts[2]);

          if(this.transcriptPositionInfo[dpo] === undefined){
            this.transcriptPositionInfo[dpo] = [];
          }
          this.transcriptPositionInfo[dpo].push([+ts[1], +ts[2], ts[3]]);


          const tInfo = {
            transcriptId: this.transcriptId(ts),
            transcriptName: ts[3],
            txStart: +ts[1],
            txEnd: +ts[2],
            strand: ts[5],
            displayOrder: dpo,
            //displayOrder: displayOrder,
          };
          this.transcriptInfo[tInfo.transcriptId] = tInfo;
          displayOrder += 1;
          // this.numTranscriptRows = Math.max(
          //   this.numTranscriptRows,
          //   displayOrder
          // );
        });

      this.numTranscriptRows = Object.keys(this.transcriptPositionInfo).length;

      //console.log(visibleTranscripts);
      //console.log(this.transcriptPositionInfo);

      //console.log(this.transcriptInfo);
    }

    calculateTranscriptRowNumber(transcriptPositionInfo,txStart,txEnd){

      const numRows = Object.keys(transcriptPositionInfo).length;

      // if(numRows === 0){
      //   return 0;
      // }
      //console.log(transcriptPositionInfo);
      //console.log(numRows);

      for(let row = 0; row < numRows; row++){
        
        let spaceAvailableOnCurrentRow = true
        transcriptPositionInfo[row]
        .forEach((ts) => {
          const currentTsStart = ts[0];
          const currentTsEnd = ts[1];

          //console.log(row, currentTsStart, currentTsEnd);
          if(
            (currentTsStart <= txStart && txStart <= currentTsEnd) ||
            (currentTsStart <= txEnd && txEnd <= currentTsEnd)
            ){
              spaceAvailableOnCurrentRow = false;
            }

        });

        if(spaceAvailableOnCurrentRow){
          return row;
        }
      }

      // If we are here, there are now available space on the existing rows.
      // Add a new one.
      return numRows;

      // transcriptPositionInfo
      //   .sort(function (a, b) {
      //     return +a[0] - b[0];
      //   })
      //   .forEach((ts) => {
      //     const currentRow = ts[0];
      //     const currentTsStart = ts[1];
      //     const currentTsEnd = ts[2];

      //     //if(txStart > )

      //   });
    }

    /*
     * Redraw the track because the options
     * changed
     */
    rerender(options, force) {
      const strOptions = JSON.stringify(options);
      if (!force && strOptions === this.prevOptions) return;
      console.log("rerender");
      //super.rerender(options, force);

      this.initOptions();

      this.prevOptions = strOptions;

      renderToggleBtn(this);

      this.updateTranscriptInfo();

      // Adjusting the track height leads to a full rerender.
      // No need to rerender again
      if(this.adjustTrackHeight()) return;

      this.visibleAndFetchedTiles().forEach((tile) => {
        this.renderTile(tile);
      });
    }

    drawTile() {}

    transcriptId(geneInfo) {
      return `${geneInfo[7]}_${geneInfo[0]}_${geneInfo[1]}_${geneInfo[2]}`;
    }

    renderTile(tile) {
      if (!tile.initialized) return;

      console.log("renderT");

      tile.allRects = [];
      // store the scale at while the tile was drawn at so that
      // we only resize it when redrawing
      tile.drawnAtScale = this._xScale.copy();
      tile.rectGraphics.removeChildren();
      tile.rectGraphics.clear();
      tile.textBgGraphics.clear();

      if(this.areTranscriptsHidden) return;


      const GENE_ALPHA = 1;

      // let plusFillerRects = tile.tileData.filter(
      //   td => td.type === 'filler' && td.strand === '+'
      // );
      // let minusFillerRects = tile.tileData.filter(
      //   td => td.type === 'filler' && td.strand === '-'
      // );

      const plusGenes = tile.tileData.filter(
        (td) =>
          td.type !== "filler" && (td.strand === "+" || td.fields[5] === "+")
      );
      const minusGenes = tile.tileData.filter(
        (td) =>
          td.type !== "filler" && (td.strand === "-" || td.fields[5] === "-")
      );

      //flagOverlappingFillers(plusGenes, plusFillerRects);
      //flagOverlappingFillers(minusGenes, minusFillerRects);

      // remove the fillers that are contained within a gene
      // plusFillerRects = plusFillerRects.filter(x => !x.hide);
      // minusFillerRects = minusFillerRects.filter(x => !x.hide);
      const yMiddle = this.transcriptHeight + this.transcriptSpacing; //this.dimensions[1] / 2;

      // const fillerGeneSpacing = (this.options.fillerHeight - this.transcriptHeight) / 2;
      const strandCenterY =this.transcriptHeight / 2 + this.transcriptSpacing / 2;

      const plusRenderContext = [
        this,
        tile,
        tile.rectGraphics,
        this._xScale,
        this.colors["+"],
        GENE_ALPHA,
        strandCenterY,
        this.transcriptHeight,
        this.transcriptSpacing,
      ];
      const minusRenderContext = [
        this,
        tile,
        tile.rectGraphics,
        this._xScale,
        this.colors["-"],
        GENE_ALPHA,
        strandCenterY,
        this.transcriptHeight,
        this.transcriptSpacing,
      ];

      //renderRects(plusFillerRects, ...plusRenderContext);
      //renderRects(minusFillerRects, ...minusRenderContext);

      // plusRenderContext[5] = GENE_ALPHA;
      // minusRenderContext[5] = GENE_ALPHA;

      renderGenes(plusGenes, ...plusRenderContext);
      renderGenes(minusGenes, ...minusRenderContext);

      renderMask(this, tile);

      trackUtils.stretchRects(this, [
        (x) => x.rectGraphics,
        (x) => x.rectMaskGraphics,
      ]);

      for (const text of Object.values(tile.texts)) {
        text.style = {
          fontSize: `${this.fontSize}px`,
          fontFamily: this.options.fontFamily,
          fill: this.colors["labelFont"],
        };
      }
      //this.draw();
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

      return zoomLevel;
    }

    draw() {
      //super.draw();

      this.allTexts = [];
      this.allBoxes = [];
      const allTiles = [];

      this.geneAreaHeight = this.transcriptHeight;
      const fontSizeHalf = this.fontSize / 2;

      trackUtils.stretchRects(this, [
        (x) => x.rectGraphics,
        (x) => x.rectMaskGraphics,
      ]);

      Object.values(this.fetchedTiles)
        // tile hasn't been drawn properly because we likely got some
        // bogus data from the server
        .filter((tile) => tile.drawnAtScale)
        .forEach((tile) => {
          tile.textBgGraphics.clear();
          tile.textBgGraphics.beginFill(
            typeof this.options.labelBackgroundColor !== "undefined"
              ? this.colors["labelBackground"]
              : WHITE_HEX
          );

          // move the texts
          //const parentInFetched = this.parentInFetched(tile);

          if (!tile.initialized) return;
          //console.log('---');
          tile.tileData.forEach((td) => {
            // tile probably hasn't been initialized yet
            if (!tile.texts) return;

            const geneInfo = td.fields;
            //console.log(geneInfo);
            const geneName = geneInfo[3];
            const geneId = this.transcriptId(geneInfo);

            if(this.transcriptInfo[geneId] === undefined) return;

            const text = tile.texts[geneId];

            if (!text) return;

            if(this.areTranscriptsHidden){
              text.visible = false;
              return;
            }

            if (!tile.textWidths[geneId]) {
              // if we haven't measured the text's width in renderTile, do it now
              // this can occur if the same gene is in more than one tile, so its
              // dimensions are measured for the first tile and not for the second
              const textWidth = text.getBounds().width;
              const textHeight = text.getBounds().height;

              tile.textHeights[geneId] = textHeight;
              tile.textWidths[geneId] = textWidth;
            }

            const TEXT_MARGIN = 3;
            const chrOffset = +td.chrOffset;
            const txStart = +geneInfo[1] + chrOffset;
            const txEnd = +geneInfo[2] + chrOffset;
            // //const txMiddle = (txStart + txEnd) / 2;
            // const txMiddle = Math.max(
            //   this.xScale().domain()[0],
            //   txStart - tile.textWidths[geneId]/2);
            //console.log(txStart, txMiddle);
            let textYMiddleOffset =
              this.transcriptInfo[geneId].displayOrder *
              (this.geneAreaHeight + this.transcriptSpacing);

            if (this.options.showToggleTranscriptsButton) {
              textYMiddleOffset += this.toggleButtonHeight;
            }

            //console.log(txStart,txEnd,txMiddle, chrOffset, tile.textHeights[geneId]);
            let textYMiddle =
              this.geneAreaHeight / 2 +
              this.transcriptSpacing / 2 +
              textYMiddleOffset;

            //const fontRectPadding = (this.geneAreaHeight - this.fontSize) / 2;


            // take care of label positioning at start or end of transcripts
            text.position.x = Math.max(
              this._xScale(this.xScale().domain()[0]) + TEXT_MARGIN,
              this._xScale(txStart) - tile.textWidths[geneId] - 2 * TEXT_MARGIN
            );

            const marginRight = this.transcriptInfo[geneId].strand === "+"
            ? tile.textWidths[geneId] + this.geneAreaHeight / 2 + 2 * TEXT_MARGIN
            : tile.textWidths[geneId] + TEXT_MARGIN

            text.position.x = Math.min(
              text.position.x,
              this._xScale(txEnd) - marginRight
            );
 
            text.position.y = textYMiddle;

            // Determine if the current text should be hidden
            let showText = true;
            //console.log(this.transcriptInfo[geneId]);
            const dpo = this.transcriptInfo[geneId].displayOrder
            //console.log(this.transcriptPositionInfo[dpo]);

            this.transcriptPositionInfo[dpo]
              .filter(ts => {
                // Check the ones that are left of the current transcript
                return ts[1] < this.transcriptInfo[geneId].txStart;
              }).forEach(ts => {
                const endOfTranscript = this._xScale(ts[1]);
                //console.log(text.position.x,endOfTranscript, tile.textWidths[geneId]);
                if(endOfTranscript > text.position.x - 4 * TEXT_MARGIN){
                  showText = false;
                }
              });


            //if (!parentInFetched) {
            if(showText){
              text.visible = true;

              this.allBoxes.push([
                text.position.x - TEXT_MARGIN,
                textYMiddle,
                tile.textWidths[geneId] + 2 * TEXT_MARGIN,
                this.transcriptHeight,
                geneName,
              ]);

              this.allTexts.push({
                importance: +geneInfo[4],
                text,
                caption: geneName,
                strand: geneInfo[5],
              });

              allTiles.push(tile.textBgGraphics);
            } else {
              text.visible = false;
            }
          });
        });

      //this.hideOverlaps(this.allBoxes, this.allTexts);
      this.renderTextBg(this.allBoxes, this.allTexts, allTiles);
    }

    renderTextBg(allBoxes, allTexts, allTiles) {
      allTexts.forEach((text, i) => {
        if (text.text.visible && allBoxes[i] && allTiles[i]) {
          const [minX, minY, width, height] = allBoxes[i];

          allTiles[i].drawRect(minX, minY - height / 2, width, height);
        }
      });
    }

    hideTexts(allTexts) {
      allTexts.forEach((text, i) => {
        text.visible = false;
      });
    }

    // hideOverlaps(allBoxes, allTexts) {
    //   boxIntersect(allBoxes, (i, j) => {
    //     if (allTexts[i].importance > allTexts[j].importance) {
    //       allTexts[j].text.visible = false;
    //     } else {
    //       allTexts[i].text.visible = false;
    //     }
    //   });
    // }

    setPosition(newPosition) {
      super.setPosition(newPosition);

      [this.pMain.position.x, this.pMain.position.y] = this.position;
    }

    setDimensions(newDimensions) {

      this.updateTranscriptInfo();

      // This will rerender all tiles.
      super.setDimensions(newDimensions);
      //console.log("setDimensions", this.dimensions);
    }

    zoomed(newXScale, newYScale) {
      this.xScale(newXScale);
      this.yScale(newYScale);

      this.refreshTiles();

      this.draw();
    }

    getMouseOverHtml(trackX, trackY) {
      if (!this.tilesetInfo) {
        return "";
      }

      const point = [trackX, trackY];

      for (const tile of this.visibleAndFetchedTiles()) {
        for (let i = 0; i < tile.allRects.length; i++) {
          // console.log('tile.allRects:', tile.allRects);
          // copy the visible rects array
          if (tile.allRects[i][2].type === "filler") {
            continue;
          }
          const rect = tile.allRects[i][0].slice(0);
          // console.log('rect:', rect);

          const newArr = [];
          while (rect.length) {
            const newPoint = rect.splice(0, 2);
            newPoint[0] =
              newPoint[0] * tile.rectGraphics.scale.x +
              tile.rectGraphics.position.x;
            newPoint[1] =
              newPoint[1] * tile.rectGraphics.scale.y +
              tile.rectGraphics.position.y;

            newArr.push(newPoint);
          }

          const pc = classifyPoint(newArr, point);

          if (pc === -1) {
            // console.log('ar:', tile.allRects[i]);
            const gene = tile.allRects[i][2];

            return `
              <div>
                <p><b>${gene.fields[3]}</b></p>
                <p>${gene.fields[0]}:${gene.fields[1]}-${gene.fields[2]} Strand: ${gene.strand}</p>
              </div>
            `;
          }
        }
      }

      return "";
    }
  }
  return new TranscritpsTrackClass(...args);
};

const icon =
  '<svg width="20" height="20" xmlns="http://www.w3.org/2000/svg"><path fill="#fff" d="M-1-1h22v22H-1z"/><g><path stroke="#007fff" stroke-width="1.5" fill="#007fff" d="M-.667-.091h5v20.167h-5z"/><path stroke-width="1.5" stroke="#e8e500" fill="#e8e500" d="M5.667.242h5v20.167h-5z"/><path stroke-width="1.5" stroke="#ff0038" fill="#ff0038" d="M15.833.076h5v20.167h-5z"/><path stroke="green" stroke-width="1.5" fill="green" d="M10.833-.258H14.5v20.167h-3.667z"/></g></svg>';

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
    "showToggleTranscriptsButton",
    "sequenceData"
  ],
  defaultOptions: {
    fontSize: 9,
    fontFamily: "Helvetica",
    transcriptSpacing: 2,
    transcriptHeight: 11,
    maxTexts: 100,
    plusStrandColor: "#bdbfff",
    minusStrandColor: "red",
    utrColor: "#C0EAAF",
    labelBackgroundColor: "#ffffff",
    labelFontColor: "#333333",
    showToggleTranscriptsButton: true
  },
};

export default TranscritpsTrack;
