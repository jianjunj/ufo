/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_CRTMIXEDSKY_H_
#define TEST_UFO_CRTMIXEDSKY_H_

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/base/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/missingValues.h"
#include "oops/util/Range.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsTraits.h"

namespace ufo {
namespace test {

void testMixedSky() {
  typedef ::test::ObsTestsFixture<ObsTraits> Fixture;
  typedef oops::ObsOperator<ObsTraits> Operator;
  typedef oops::ObsVector<ObsTraits> Vector;
  const auto conf = ::test::TestEnvironment::config().getSubConfigurations("observations")[0];
  const auto & space = Fixture::obspace()[0];

  Operator mixed(space, conf.getSubConfiguration("obs operator"));
  Operator cloudy(space, conf.getSubConfiguration("cloudy obs operator"));
  Operator clear(space, conf.getSubConfiguration("clear obs operator"));
  oops::ObsAuxControl<ObsTraits> biasControl(space, conf.getSubConfiguration("obs bias"));
  oops::Variables vars = mixed.requiredVars();
  vars += cloudy.requiredVars();
  vars += clear.requiredVars();
  vars += biasControl.requiredVars();
  oops::GeoVaLs<ObsTraits> geovals(conf.getSubConfiguration("geovals"), space, vars);
  mixed.computeReducedVars(biasControl.requiredVars(), geovals);
  oops::ObsDataVector<ObsTraits, int> qc(space, space.obsvariables(), std::string());
  Vector mixedValues(space), cloudyValues(space), clearValues(space), bias(space);
  bias.zero();
  oops::ObsVariables diagvars;
  diagvars += biasControl.requiredHdiagnostics();
  oops::ObsDiagnostics<ObsTraits> diags(space, mixed.locations(), diagvars);
  mixed.simulateObs(geovals, mixedValues, biasControl, qc, bias, diags);
  cloudy.simulateObs(geovals, cloudyValues, biasControl, qc, bias, diags);
  clear.simulateObs(geovals, clearValues, biasControl, qc, bias, diags);

  const auto channels = space.assimvariables().channels();
  const auto & mixedVec = mixedValues.obsvector();
  const auto & cloudyVec = cloudyValues.obsvector();
  const auto & clearVec = clearValues.obsvector();
  const double missing = util::missingValue<double>();
  size_t checked = 0;
  size_t changed = 0;
  for (size_t loc = 0; loc < mixedVec.nlocs(); ++loc) {
    for (size_t ch = 0; ch < channels.size(); ++ch) {
      const size_t index = loc * channels.size() + ch;
      const bool clearChannel = channels[ch] == 1 || channels[ch] == 4 || channels[ch] >= 14;
      const double expected = clearChannel ? clearVec[index] : cloudyVec[index];
      if (expected == missing) {
        EXPECT(mixedVec[index] == missing);
      } else {
        EXPECT(std::abs(mixedVec[index] - expected) < 1.0e-7);
        ++checked;
        if (clearChannel && cloudyVec[index] != missing &&
            std::abs(clearVec[index] - cloudyVec[index]) > 1.0e-5) ++changed;
      }
    }
  }
  EXPECT(checked > 0);
  EXPECT(changed > 0);

  std::vector<int> zeroCloud(mixedVec.nlocs(), 0);
  for (size_t loc = 0; loc < zeroCloud.size(); loc += 2) zeroCloud[loc] = 1;
  space.obsspace().put_db("MetaData", "zeroCloudInCRTM", zeroCloud);
  Vector maskedValues(space);
  mixed.simulateObs(geovals, maskedValues, biasControl, qc, bias, diags);
  const auto & maskedVec = maskedValues.obsvector();
  for (size_t loc = 0; loc < mixedVec.nlocs(); ++loc) {
    for (size_t ch = 0; ch < channels.size(); ++ch) {
      const size_t index = loc * channels.size() + ch;
      const double expected = zeroCloud[loc] ? clearVec[index] : mixedVec[index];
      if (expected == missing) {
        EXPECT(maskedVec[index] == missing);
      } else {
        EXPECT(std::abs(maskedVec[index] - expected) < 1.0e-6);
      }
    }
  }
}

void testMixedSkyTLADMask() {
  typedef ::test::ObsTestsFixture<ObsTraits> Fixture;
  const auto conf = ::test::TestEnvironment::config().getSubConfigurations("observations")[0];
  const auto & space = Fixture::obspace()[0];
  const auto channels = space.assimvariables().channels();
  std::vector<int> zeroCloud(space.obsspace().nlocs(), 0);
  for (size_t loc = 0; loc < zeroCloud.size(); loc += 2) zeroCloud[loc] = 1;
  space.obsspace().put_db("MetaData", "zeroCloudInCRTM", zeroCloud);

  ufo::LinearObsOperator linear(space.obsspace(),
                                conf.getSubConfiguration("obs operator"));
  ufo::GeoVaLs trajectory(conf.getSubConfiguration("geovals"), space.obsspace(),
                          linear.requiredVars());
  ufo::ObsBias linearBias(space.obsspace(),
                          conf.getSubConfiguration("obs bias"));
  ufo::ObsBiasIncrement biasIncrement(space.obsspace(),
                                      conf.getSubConfiguration("obs bias"));
  ioda::ObsDataVector<int> qc(space.obsspace(), space.obsspace().obsvariables());
  linear.setTrajectory(trajectory, linearBias, qc);

  ufo::GeoVaLs cloudWaterIncrement(trajectory);
  cloudWaterIncrement.zero();
  const oops::Variable cloudWater(
      "mass_content_of_cloud_liquid_water_in_atmosphere_layer");
  std::vector<double> cloudWaterProfile(cloudWaterIncrement.nlevs(cloudWater),
                                        1.0e-4);
  for (size_t profile = 0; profile < cloudWaterIncrement.nprofiles(cloudWater);
       ++profile) {
    cloudWaterIncrement.putProfile(cloudWaterProfile, cloudWater, profile);
  }

  ioda::ObsVector tlValues(space.obsspace());
  tlValues.zero();
  linear.simulateObsTL(cloudWaterIncrement, tlValues, biasIncrement);

  ioda::ObsVector adjointValues(space.obsspace());
  adjointValues.zero();
  const auto channelTwo = std::find(channels.begin(), channels.end(), 2);
  ASSERT(channelTwo != channels.end());
  const size_t channelTwoIndex = std::distance(channels.begin(), channelTwo);
  for (size_t loc = 0; loc < zeroCloud.size(); ++loc) {
    adjointValues[loc * channels.size() + channelTwoIndex] = 1.0;
  }
  ufo::GeoVaLs cloudWaterAdjoint(trajectory);
  cloudWaterAdjoint.zero();
  linear.simulateObsAD(cloudWaterAdjoint, adjointValues, biasIncrement);

  std::vector<util::Range<size_t>> profilesByLocation;
  cloudWaterAdjoint.getProfileIndicesGroupedByLocation(cloudWater,
                                                       profilesByLocation);
  double unmaskedTL = 0.0;
  double unmaskedAD = 0.0;
  for (size_t loc = 0; loc < zeroCloud.size(); ++loc) {
    const size_t channelIndex = loc * channels.size() + channelTwoIndex;
    if (zeroCloud[loc]) {
      for (size_t ch = 0; ch < channels.size(); ++ch) {
        EXPECT(std::abs(tlValues[loc * channels.size() + ch]) < 1.0e-12);
      }
    } else {
      unmaskedTL += std::abs(tlValues[channelIndex]);
    }

    double adMagnitude = 0.0;
    for (size_t profile = profilesByLocation[loc].begin;
         profile < profilesByLocation[loc].end; ++profile) {
      std::vector<double> profileValues;
      cloudWaterAdjoint.getProfile(profileValues, cloudWater, profile);
      for (const double value : profileValues) adMagnitude += std::abs(value);
    }
    if (zeroCloud[loc]) {
      EXPECT(adMagnitude < 1.0e-12);
    } else {
      unmaskedAD += adMagnitude;
    }
  }
  EXPECT(unmaskedTL > 1.0e-12);
  EXPECT(unmaskedAD > 1.0e-12);
}

class CRTMixedSky : public oops::Test {
  typedef ::test::ObsTestsFixture<ObsTraits> Fixture;
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "ufo::test::CRTMixedSky";}
  void register_tests() const override {
    eckit::testing::specification().emplace_back(CASE("ufo/CRTMixedSky")
      { testMixedSky(); });
    eckit::testing::specification().emplace_back(CASE("ufo/CRTMixedSky/TLADMask")
      { testMixedSkyTLADMask(); });
  }
  void clear() const override {Fixture::reset();}
};

}  // namespace test
}  // namespace ufo
#endif  // TEST_UFO_CRTMIXEDSKY_H_
