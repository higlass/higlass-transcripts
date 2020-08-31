import { expect } from "chai";
import register from "higlass-register";

import FetchMockHelper from "./utils/FetchMockHelper";

import { HiGlassComponent, getTrackObjectFromHGC } from "higlass";

import {
  waitForDataLoaded,
  mountHGComponent,
  removeHGComponent,
} from "./utils/test-helpers";

import viewConf from "./view-configs/simple-track";

import TranscriptsTrack from "../src/scripts/TranscriptsTrack";

register({
  name: "TranscriptsTrack",
  track: TranscriptsTrack,
  config: TranscriptsTrack.config,
});

describe("SVG export", () => {
  const fetchMockHelper = new FetchMockHelper("", "SVGExport");

  beforeAll(async () => {
    await fetchMockHelper.activateFetchMock();
  });

  describe("SVG export", () => {
    let hgc = null;
    let div = null;

    beforeAll((done) => {
      [div, hgc] = mountHGComponent(div, hgc, viewConf, done);
    });

    it("tests that the export works and contains the correct data", (done) => {
      hgc.instance().handleExportSVG();

      const trackObj = getTrackObjectFromHGC(
        hgc.instance(),
        viewConf.views[0].uid,
        viewConf.views[0].tracks.top[1].uid
      );

      const tile = trackObj.visibleAndFetchedTiles()[1];
      const exon = tile.allExonsForSVG[1];

      expect(exon.rect[0]).to.equal(285.4775);
      expect(exon.rect[5]).to.equal(33);
      expect(exon.color).to.equal("#C0EAAF");

      const exon2 = tile.allExonsForSVG[2];

      expect(exon2.rect[0]).to.equal(324.7475);
      expect(exon2.rect[5]).to.equal(33);
      expect(exon2.color).to.equal("#bdbfff");


      done();
    });

    afterAll(() => {
      removeHGComponent(div);
    });
  });

  afterAll(async () => {
    await fetchMockHelper.storeDataAndResetFetchMock();
  });
});
