#include "geogram_utils.h"


// Put geogram log messages in the trash where they belong
class GeoTrashCan: public GEO::LoggerClient {
protected:
    void div(const std::string& log) override { /*pybind11::print(log);*/ }
    void out(const std::string& log) override { /*pybind11::print(log);*/ }
    void warn(const std::string& log) override { /*pybind11::print(log);*/ }
    void err(const std::string& log) override { /*pybind11::print(log);*/ }
    void status(const std::string& log) override { /*pybind11::print(log);*/ }
};


// Geogram is stateful and needs to be initialized.
// These variables keep track of whether geogram is initialized in a thread-safe manner.
// I'm using p-threads so none of this will work on Windows.
bool geogram_is_initialized = false;
std::mutex geogram_init_mutex;

// Initialize geogram exactly once.
void init_geogram_only_once() {
  std::lock_guard<std::mutex> guard(geogram_init_mutex);

  if (!geogram_is_initialized) {
    GEO::initialize();

    GEO::Logger *geo_logger = GEO::Logger::instance();
    geo_logger->unregister_all_clients();
    geo_logger->register_client(new GeoTrashCan());
    geo_logger->set_pretty(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    geogram_is_initialized = true;
  }
}
