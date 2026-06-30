#include "stage1_beam_mapper.h"
#include "stage1_config.h"
#include "stage1_data_reader.h"
#include "stage1_pulse_resampler.h"
#include "stage1_range_resizer.h"
#include "stage1_validator.h"
#include "stage1_writer.h"

#include <exception>
#include <iostream>

using namespace gmti::sim_stage1;

int main(int argc, char **argv)
{
    Stage1RunOptions opt;
    Stage1NewSystemConfig new_cfg;
    std::string err;
    try {
        parseCommandLine(argc, argv, opt, new_cfg);
    } catch (const std::exception &e) {
        std::cerr << "[stage1][ERR] " << e.what() << "\n";
        return 2;
    }
    if (!ensureStage1Dirs(opt.output_dir, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 2;
    }

    Stage1OldSystemConfig old_cfg;
    if (!parseOldConfigXml(opt.old_config_path, old_cfg, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 2;
    }
    if (old_cfg.az_count == 25) {
        std::cerr << "[stage1][WARN] XML az_count=25, but stage1 old beam table is defined as 26. Using 26 for audit/generation.\n";
        old_cfg.az_count = 26;
    }

    FileAudit audit;
    if (!auditOldFiles(old_cfg, audit)) {
        writeParameterAudit(old_cfg, new_cfg, audit, opt.output_dir, err);
        std::cerr << "[stage1][ERR] file audit failed: " << audit.message << "\n";
        return 3;
    }

    std::vector<BeamMapEntry> beam_map;
    if (!buildBeamMap(new_cfg, beam_map, err) ||
        !writeBeamTables(opt.output_dir, beam_map, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 4;
    }

    std::vector<PulseMapEntry> pulse_map;
    std::vector<PulseMapEntry> pulse_map_norm;
    if (!buildPulseMap(old_cfg, new_cfg, opt.pulse_resample_mode, pulse_map, err) ||
        !writePulseMap(opt.output_dir, opt.pulse_resample_mode, pulse_map, err) ||
        !buildPulseMap(old_cfg, new_cfg, "normalized", pulse_map_norm, err) ||
        !writePulseMap(opt.output_dir, "normalized", pulse_map_norm, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 5;
    }
    if (opt.pulse_resample_mode != "physical_time") {
        std::vector<PulseMapEntry> physical;
        if (!buildPulseMap(old_cfg, new_cfg, "physical_time", physical, err) ||
            !writePulseMap(opt.output_dir, "physical_time", physical, err)) {
            std::cerr << "[stage1][ERR] " << err << "\n";
            return 5;
        }
    }

    if (!writeParameterAudit(old_cfg, new_cfg, audit, opt.output_dir, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 6;
    }

    RangeFftZeroPadResizer resizer(old_cfg.pulse_len, new_cfg.ddc_len_new);
    if (opt.validate && !runRangeResizeSelfTests(opt.output_dir, resizer, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 7;
    }

    const std::string out_file = pathJoin(pathJoin(opt.output_dir, "data"), "stage1_background_newprotocol.bin");
    GenerationStats stats;
    if (opt.generate_data) {
        if (!generateStage1Data(old_cfg, new_cfg, opt, beam_map, pulse_map,
                                resizer, out_file, stats, err)) {
            std::cerr << "[stage1][ERR] " << err << "\n";
            return 8;
        }
    }
    if (opt.write_config && !writeStage1ConfigXml(old_cfg, new_cfg, opt, out_file, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 9;
    }
    if (!writeDataStatsReport(opt.output_dir, stats, opt, new_cfg, err) ||
        !writeAlgorithmIntegrationPlaceholder(opt.output_dir, opt, new_cfg, stats, err)) {
        std::cerr << "[stage1][ERR] " << err << "\n";
        return 10;
    }

    std::cout << "[stage1] done. output_dir=" << opt.output_dir << "\n";
    return 0;
}

