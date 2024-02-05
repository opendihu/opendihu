#pragma once

#include <Python.h> // has to be the first included header

#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <tuple>
#include <map>

//! forward declarations
namespace Data {
template <typename FunctionSpaceType, int nComponents1, int nComponents2>
class SlotConnectorData;
}
class SlotsConnection;

/** Class that collects information about all nested solvers and produces a
 * diagram as txt file.
 *
 */
class SolverStructureVisualizer {
public:
  /** Representation of one solver in the diagram
   */
  struct solver_t {
    std::string name;        //< name (type) of the solver
    std::string description; //< additional string that will be included, e.g.
                             // for type of subsolver
    bool hasInternalConnectionToFirstNestedSolver; //< if the solver has an
                                                   // internal connection of all
                                                   // connector slots of its
                                                   // first subsolver. This is
                                                   // the case e.g. for Coupling
                                                   // and StrangSplitting.
    bool
        hasInternalConnectionToSecondNestedSolver; //< if the solver has an
                                                   // internal connection of all
                                                   // connector slots of its
                                                   // second subsolver. This is
                                                   // the case e.g. for Coupling
                                                   // and StrangSplitting.

    /** Representation of a single output slot
     */
    struct OutputSlot {
      int variableNo; //< either 0 or 1, if the slot is internally stored in
                      // variable1 or variable2 of SlotConnectorData
      std::string fieldVariableName; //< the name of the field variable that
                                     // will be written in the diagram
      std::string componentName;     //< the name of the component that will be
                                     // written in the diagram
      int nComponents; //< number of components the field variable has in order
                       // to distinguish scalar field variables
      std::string meshDescription; //< string information of the mesh, from
                                   // FunctionSpace::getDescription()
      std::string slotName;        //< slot name if given
    };

    std::vector<OutputSlot> outputSlots;

    std::shared_ptr<SlotsConnection>
        slotsConnection; //< pointer to the actual slotsConnection object of the
                         // operator splitting

    /** connection between two output slots of two children
     */
    struct SlotsConnectionRepresentation {
      int fromSlot;
      int toSlot;
      enum slot_connection_t {
        ab,
        ba,
        bidirectionalCopy,
        bidirectionalReuse,
        internalBeforeComputation,
        internalAfterComputation
      } type; //< ab=term1 -> term2, ba=term2 -> term1, bidirectional=shared
              // between term1 and term2,
              // internalBeforeComputation=MapDofs.beforeComputation,
              // internalAfterComputation=MapDofs.afterComputation
      bool involvesMapping;
    };

    std::vector<SlotsConnectionRepresentation>
        slotsConnections; //< connections between output slots
    std::vector<SlotsConnectionRepresentation>
        mappingsWithinSolver; //< "connections" within the same solver, this is
                              // used for MapDofs

    std::vector<std::shared_ptr<solver_t>>
        children; //< the nested solvers inside the current solver
    std::shared_ptr<solver_t>
        parent; //< pointer to the parent of the current nested solver
  };

  //! constructor
  SolverStructureVisualizer();

  //! add a solver to the diagram
  void addSolver(std::string name,
                 bool hasInternalConnectionToFirstNestedSolver = false,
                 bool hasInternalConnectionToSecondNestedSolver = false);

  //! add a description for the current that will be included in the
  //! visualization, to be called after addSolver
  void setSolverDescription(std::string description);

  //! indicate that all further calls to addSolver will be children of the
  //! current solver
  void beginChild(std::string description = "");

  //! indicate the end of the current child
  void endChild();

  //! add the output connection information between two children to the current
  //! solver
  void addSlotsConnection(std::shared_ptr<SlotsConnection> slotsConnection);

  //! set the slot connector data
  template <typename FunctionSpaceType, int nComponents1, int nComponents2>
  void setSlotConnectorData(std::shared_ptr<Data::SlotConnectorData<
                                FunctionSpaceType, nComponents1, nComponents2>>
                                slotConnectorData,
                            bool isFromTuple = false);

  //! set the slot connector data
  template <typename T>
  void setSlotConnectorData(std::shared_ptr<std::vector<T>> slotConnectorData,
                            bool isFromTuple = false);

  //! set the slot connector data
  template <typename SlotConnectorData1, typename SlotConnectorData2>
  void setSlotConnectorData(
      std::shared_ptr<std::tuple<SlotConnectorData1, SlotConnectorData2>>
          slotConnectorData,
      bool isFromTuple = false);

  //! add connections between slots that occur within the same solver, this is
  //! used by MapDofs
  void addSlotMapping(
      int slotNoFrom, int slotNoTo,
      solver_t::SlotsConnectionRepresentation::slot_connection_t internalType);

  //! beginChild, addSolver etc. as well as writeDiagramFile will have no
  //! effect, this is needed for MultipleInstances such that not all instances
  //! appear as separate solvers but only the first one
  void disable();

  //! beginChild and addSolver and writeDiagramFile etc. will work again
  void enable();

  //! produce the resulting file, only if solverAddingEnabled_
  void writeDiagramFile(std::string filename);

  //! get the diagram as string
  std::string getDiagram();

protected:
  /** a nested class that recursively traverses the tree of nested solvers and
   * generates the diagram
   */
  struct DiagramGenerator {
    //! initialize the diagram generator, afterwards, generateDiagram can be
    //! called
    void initialize(std::shared_ptr<solver_t> solverRoot);

    //! generate the diagram
    std::string generateFinalDiagram();

  private:
    //! print the nested solver structure to result
    //! @param externalConnectionLines <fromLineNo, toLineNo, linePosition,
    //! type>
    void generateDiagramRecursion(
        std::stringstream &result,
        std::vector<std::vector<int>> &internalConnectionLinesCurrentSolver,
        std::vector<int> &slotLineNos, int depth, bool isFirstChildSolver,
        bool isLastChildSolver);

    //! generate a diagram without the vertical connection lines, neither
    //! internal nor external lines
    std::string generateDiagramWithoutConnectionLines();

    //! internal variables that are used by the enclosing class to add
    //! connection lines to the diagram
    std::shared_ptr<solver_t>
        currentSolver_; //< a pointer to the current solver during the iteration

    /** Type to store information about an external connection line between two
     * slots.
     */
    struct ExternalConnectionLine {
      int lineNoFrom; //< row where the line starts
      int lineNoTo;   //< row where the line ends
      int lineColumn; //< lineColumn is the horizontal position of the vertical
                      // data connection line
      int slotNoFrom; //< the slot no where the line starts, this is needed to
                      // know the number of digits of the slot no and where to
                      // put the '<' sign
      int slotNoTo; //< the slot no where the line ends, this is needed to know
                    // the number of digits of the slot no and where to put the
                    //'<' sign
      SolverStructureVisualizer::solver_t::SlotsConnectionRepresentation::
          slot_connection_t
              lineType;     //< type of the line if it is copy or reuse
      bool involvesMapping; //< if the line is a mapping
    };
    std::vector<ExternalConnectionLine>
        externalConnectionLines_; //< connection lines contains the following
                                  // information: <lineNoFrom, lineNoTo,
                                  // lineColumn, lineType>
    std::vector<std::vector<int>>
        internalConnectionLinesAll_; //< list of line nos that are connected by
                                     // internal connections
    std::vector<int>
        slotLineNosAll_; //< all slots in the diagram with their line no
    std::map<int, std::string>
        meshDescriptionsOfSlots_; //< for every line no of a slot, the mesh name
                                  // of that slot, used to determine if a
                                  // mapping takes place
    std::vector<std::string> referencedMeshNames_; //< all occuring mesh names
  };

  //! in the currentSolver_ fill slotsConnections vector from slotsConnection
  static void parseSlotsConnection(std::shared_ptr<solver_t> currentSolver);

  std::shared_ptr<solver_t> solverRoot_; //< the whole nested solver structure
  std::shared_ptr<solver_t>
      currentSolver_; //< a pointer to the current solver for which call to
                      // addSolver sets the name and data

  bool enabled_;      //< if addSolver has an effect
  int nDisableCalls_; //< how often disable() has been called in sequence
};

#include "control/diagnostic_tool/solver_structure_visualizer.tpp"
