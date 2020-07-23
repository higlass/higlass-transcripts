import { expect } from "chai";
import register from "higlass-register";

import FetchMockHelper from "./utils/FetchMockHelper";

import { HiGlassComponent, getTrackObjectFromHGC } from "higlass";

import {
  waitForDataLoaded,
  mountHGComponent,
  removeHGComponent,
} from "./utils/test-helpers";

import viewConf from "./view-configs/simple-track-aa";

import TranscriptsTrack from "../src/scripts/TranscriptsTrack";

register({
  name: "TranscriptsTrack",
  track: TranscriptsTrack,
  config: TranscriptsTrack.config,
});

describe("Amino acid", () => {
  const fetchMockHelper = new FetchMockHelper("", "AminoAcidTests");

  beforeAll(async () => {
    await fetchMockHelper.activateFetchMock();
  });

  describe("Amino acid tests", () => {
    let hgc = null;
    let div = null;

    beforeAll((done) => {
      [div, hgc] = mountHGComponent(div, hgc, viewConf, done);
    });

    it("tests that the export works and contains the correct data", (done) => {

      const trackObj = getTrackObjectFromHGC(
        hgc.instance(),
        viewConf.views[0].uid,
        viewConf.views[0].tracks.top[0].uid
      );
      const tile = trackObj.visibleAndFetchedTiles()[0];

      setTimeout(() => {  
        
        
        console.log(tile.aaInfo.exonOffsets["ENST00000378404.3_chr1_3021483_3022903"]);
        console.log(tile.aaInfo.aminoAcids);

        done();
       }, 2000);





     


      
    });

    afterAll(() => {
      removeHGComponent(div);
    });
  });

  afterAll(async () => {
    await fetchMockHelper.storeDataAndResetFetchMock();
  });
});
