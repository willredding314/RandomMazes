
import java.awt.*;
import java.util.*;

import tester.*;
import javalib.impworld.*;
import javalib.worldimages.*;

import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;

// Represents all constants used throughout the maze construction and searching
final class Constants {
  static final Color CELL_SEEN_COLOR = Color.CYAN;
  static final Color CELL_UNSEEN_COLOR = Color.LIGHT_GRAY;
  static final Color CELL_PATH_COLOR = Color.BLUE;
  static final Color CELL_START_COLOR = Color.RED;
  static final Color CELL_FINISH_COLOR = Color.MAGENTA;
  static final Color CELL_HEAD_COLOR = Color.YELLOW;
  static final Color CELL_FUTUREHEAD_COLOR = Color.ORANGE;

  static final Color CELL_CLOSE_COLOR = Color.BLUE;
  static final Color CELL_FAR_COLOR = Color.RED;

  static final Color WALL_COLOR = Color.BLACK;
  static final Color INTERPOLATE_WALL_COLOR = Color.BLACK;

  static final int WINDOW_WIDTH = 800;
  static final int WINDOW_HEIGHT = 800;

  // Prevent construction. Should only
  // access the static fields of the class
  private Constants() {
  }
}

// Represents a maze as a grid of cells, with edges
// connecting every side of the cell to another cell
class Maze {
  private final int cellsPerRow;
  private final int cellsPerColumn;

  // Whether or not the cells in this maze
  // are hexagonal
  private final boolean isHexagonalMaze;

  // All cells in the grid, not including DummyCells representing the boundary
  private final ArrayList<ArrayList<ICell>> cells;

  // All edges connecting Cells (does not include edges connecting to DummyCells).
  // These edges alias the same edges referred to by each ICell in the list above
  private ArrayList<CellEdge> edgesInGrid;

  // an object that creates new random weights for edges
  private final EdgeWeightGen weightMaker;

  // All edges in the minimum spanning tree after implementing Kruskal's
  // algorithm to connect all cells with open pathways
  private ArrayList<CellEdge> edgesInMST;

  // constructs a maze with a given width and height, a seed to generate the maze,
  // and two booleans to determine if the maze will be biased towards rows or
  // columns
  Maze(int cellsPerRow, int cellsPerColumn, int seed, boolean biasedH, boolean biasedV,
      boolean isHexagonalMaze) {

    // Ensure that there is at least 1 cell in each row and column
    if (cellsPerRow <= 0) {
      throw new IllegalArgumentException(
          "A maze must have at least one cell per row. Given " + Integer.toString(cellsPerRow));
    }

    if (cellsPerColumn <= 0) {
      throw new IllegalArgumentException("A maze must have at least one cell per column. Given "
          + Integer.toString(cellsPerColumn));
    }

    this.cellsPerColumn = cellsPerColumn;
    this.cellsPerRow = cellsPerRow;
    this.isHexagonalMaze = isHexagonalMaze;

    // determine what edge bias we are using, and assign weights to edges
    // accordingly
    // Note: for all edges, we assign the seeds for horizontal and vertical weights
    // as _seed_ and _seed_ + 1. The "+ 1" is arbitrary, as we only need to ensure
    // that the two seeds are different, so any change to the value would work.
    if (!biasedH && !biasedV) {

      // with no bias, we have the same maximum weight for both vertical
      // and horizontal edges
      this.weightMaker = new EdgeWeightGen(1000, 1000, seed, seed + 1);

    } else if (!biasedH && biasedV) {

      // with vertical bias, we drop the maximum vertical weight
      this.weightMaker = new EdgeWeightGen(1000, 10, seed, seed + 1);

    } else if (biasedH && !biasedV) {

      // with horizontal bias, we drop the maximum horizontal
      this.weightMaker = new EdgeWeightGen(10, 1000, seed, seed + 1);
    } else {

      // if we claim to be both vertically and horizontally biased, throw an
      // exception
      throw new IllegalArgumentException(
          "A maze cannot be both horizontally " + "and vertically biased at the same time");
    }

    if (isHexagonalMaze) {

      // ACCUMULATOR: Keeps track of all of the cells built for the maze
      ArrayList<ArrayList<ICell>> hiddenCells = new ArrayList<ArrayList<ICell>>();

      // ACCUMULATOR: Keeps track of all of the cells built for the maze
      // while keeping their identities. The cells in this list alias
      // those in _hiddenCells_ above in the same order, and are thereby
      // in a one-to-one correspondence
      ArrayList<ArrayList<HexCell>> cells = new ArrayList<ArrayList<HexCell>>();

      // Create and fill the grid of cells, initially disconnected
      // from each other but connected to 4 of their own `DummyCell`s
      for (int i = 0; i < cellsPerColumn; i += 1) {
        ArrayList<HexCell> newRow = new ArrayList<HexCell>();
        ArrayList<ICell> newHiddenRow = new ArrayList<ICell>();

        for (int n = 0; n < cellsPerRow; n += 1) {
          HexCell newCell = new HexCell();
          newRow.add(newCell);
          newHiddenRow.add(newCell);
        }
        cells.add(newRow);
        hiddenCells.add(newHiddenRow);
      }

      this.cells = hiddenCells;
      this.edgesInGrid = new ArrayList<CellEdge>();

      // With HexCells, things are slightly more complicated than before.
      // Since our Maze is a representation of hexagons turned on their sides,
      // any two adjacent rows in a maze actually represent cells shifted over
      // one cell with respect to one another. So we have something like
      //
      // ⬡⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves
      // ⬡⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves
      //
      // Thus, connections between cells occurs in the following manner:
      //
      // In this case, ⬢ connects to its bottom right and bottom left, as well as its
      // to its right neighbor
      // ⬡⬡⬡⬢⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      //
      // In this case, ⬢ connects to its bottom left and bottom right neighbors
      // ⬡⬡⬡⬡⬡⬡⬢
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      //
      // In this case, ⬢ connects to its bottom left neighbor ONLY
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬢
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      //
      // In this case, ⬢ connects to its bottom right ONLY
      // ⬢⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      //
      // In this case, ⬢ connects to its bottom left and bottom right ONLY
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬢⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      // ⬡⬡⬡⬡⬡⬡⬡
      //
      // The bottom edge only needs to connect to its right neighbors
      // since its top right and top left neighbors already connected to it

      // Connect all of the "easy" hexagons: those not on either of the right
      // or left edges.
      for (int i = 0; i < cellsPerColumn - 1; i += 1) {

        // ACCUMULATOR: Determines if we are in a "normal" row
        // or are offset with respect to the normal rows. This
        // appears diagrammatically as follows:
        // ⬡⬡⬡⬡⬡⬡⬡⬡ -> false -> row = 0 even
        // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves -> true
        // ⬡⬡⬡⬡⬡⬡⬡⬡ -> false -> row = 2 even
        // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves -> true
        //
        // We are shifted over if we are in an "odd" row starting at 0
        boolean shiftedOver = i % 2 == 1;

        // Notice we start with j = 1 here
        for (int j = 1; j < cellsPerRow - 1; j += 1) {

          int bottomLeftIndex = 0;

          // If we are in a row that is shifted with respect
          // to other rows, our bottom left is actually at the same
          // index in the subsequent row because we are mapping
          // this shifted-over structure to a grid that isn't shifted over.
          // That is, the cell at index `j` in the next row would correspond to the
          // cell that "fits" in between the `(j-1)`th cell in this row and the
          // `j`th cell in this row
          if (shiftedOver) {
            bottomLeftIndex = j;
          } else {
            bottomLeftIndex = j - 1;
          }

          HexCell currCell = cells.get(i).get(j);
          HexCell rightCell = cells.get(i).get(j + 1);
          HexCell bottomLeftCell = cells.get(i + 1).get(bottomLeftIndex);
          HexCell bottomRightCell = cells.get(i + 1).get(bottomLeftIndex + 1);

          CellEdge rightEdge = currCell.connectAtRight(rightCell,
              this.weightMaker.nextHorizontalWeight());
          CellEdge bottomLeftEdge = currCell.connectAtBottomLeft(bottomLeftCell,
              this.weightMaker.nextVerticalWeight());
          CellEdge bottomRightEdge = currCell.connectAtBottomRight(bottomRightCell,
              this.weightMaker.nextVerticalWeight());

          // Add these new edges to the edges in the maze. These edges alias
          // the actual edges connecting the cells (that the cells also reference)
          this.edgesInGrid.add(rightEdge);
          this.edgesInGrid.add(bottomLeftEdge);
          this.edgesInGrid.add(bottomRightEdge);
        }
      }

      // Connect all of the hexagons on the left-most edge
      // to their right neighbors and either both their
      // bottom-right and bottom-left neighbors or only their
      // bottom-right neighbor
      // AND
      // Connect all of the hexagons on the right-most edge
      // to either both their bottom-right and bottom-left neighbors or only their
      // bottom-left neighbor
      for (int i = 0; i < cellsPerColumn - 1; i += 1) {

        // ACCUMULATOR: Determines if we are in a "normal" row
        // or are offset with respect to the normal rows. This
        // appears diagrammatically as follows:
        // ⬡⬡⬡⬡⬡⬡⬡⬡ -> false -> row = 0 even
        // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves -> true
        // ⬡⬡⬡⬡⬡⬡⬡⬡ -> false -> row = 2 even
        // ⬡⬡⬡⬡⬡⬡⬡⬡ <-- shifted over one to fit in the grooves -> true
        //
        // We are shifted over if we are in an "odd" row starting at 0
        boolean shiftedOver = i % 2 == 1;

        // If we are shifted over with respect to other rows, then
        // we connect to both our bottom-left and bottom-right
        // neighbors
        if (shiftedOver) {

          // For cell at the left-most edge
          {
            HexCell currCell = cells.get(i).get(0);
            HexCell rightCell = cells.get(i).get(1);
            HexCell bottomLeftCell = cells.get(i + 1).get(0);
            HexCell bottomRightCell = cells.get(i + 1).get(1);

            CellEdge rightEdge = currCell.connectAtRight(rightCell,
                this.weightMaker.nextHorizontalWeight());
            CellEdge bottomLeftEdge = currCell.connectAtBottomLeft(bottomLeftCell,
                this.weightMaker.nextVerticalWeight());
            CellEdge bottomRightEdge = currCell.connectAtBottomRight(bottomRightCell,
                this.weightMaker.nextVerticalWeight());

            // Add these new edges to the edges in the maze. These edges alias
            // the actual edges connecting the cells (that the cells also reference)
            this.edgesInGrid.add(rightEdge);
            this.edgesInGrid.add(bottomLeftEdge);
            this.edgesInGrid.add(bottomRightEdge);
          }

          // For cell at the right-most edge
          {
            HexCell currCell = cells.get(i).get(cellsPerRow - 1);
            HexCell bottomLeftCell = cells.get(i + 1).get(cellsPerRow - 1);

            CellEdge bottomLeftEdge = currCell.connectAtBottomLeft(bottomLeftCell,
                this.weightMaker.nextVerticalWeight());

            // Add these new edges to the edges in the maze. These edges alias
            // the actual edges connecting the cells (that the cells also reference)
            this.edgesInGrid.add(bottomLeftEdge);
          }
        } else {
          // Otherwise, we only connect to our bottom right neighbor
          // which is the FIRST cell in the row below us because it is shifted
          // with respect to us (hence why `bottomRightCell = cells.get(i + 1).get(0))`

          // For cell at the left-most edge
          {
            HexCell currCell = cells.get(i).get(0);
            HexCell rightCell = cells.get(i).get(1);
            HexCell bottomRightCell = cells.get(i + 1).get(0);

            CellEdge rightEdge = currCell.connectAtRight(rightCell,
                this.weightMaker.nextHorizontalWeight());
            CellEdge bottomRightEdge = currCell.connectAtBottomRight(bottomRightCell,
                this.weightMaker.nextVerticalWeight());

            // Add these new edges to the edges in the maze. These edges alias
            // the actual edges connecting the cells (that the cells also reference)
            this.edgesInGrid.add(rightEdge);
            this.edgesInGrid.add(bottomRightEdge);
          }

          // For cell at the right-most edge
          {
            HexCell currCell = cells.get(i).get(cellsPerRow - 1);
            HexCell bottomRightCell = cells.get(i + 1).get(cellsPerRow - 1);
            HexCell bottomLeftCell = cells.get(i + 1).get(cellsPerRow - 2);

            CellEdge bottomLeftEdge = currCell.connectAtBottomLeft(bottomLeftCell,
                this.weightMaker.nextVerticalWeight());
            CellEdge bottomRightEdge = currCell.connectAtBottomRight(bottomRightCell,
                this.weightMaker.nextVerticalWeight());

            // Add these new edges to the edges in the maze. These edges alias
            // the actual edges connecting the cells (that the cells also reference)
            this.edgesInGrid.add(bottomLeftEdge);
            this.edgesInGrid.add(bottomRightEdge);
          }
        }
      }

      // Connect all cells in the bottom edge ONLY to their right neighbors
      for (int i = 0; i < cellsPerRow - 1; i += 1) {
        HexCell currCell = cells.get(cellsPerColumn - 1).get(i);
        HexCell rightCell = cells.get(cellsPerColumn - 1).get(i + 1);
        CellEdge rightEdge = currCell.connectAtRight(rightCell,
            this.weightMaker.nextHorizontalWeight());

        // Add the new edge to the edges in the maze. This edge aliases
        // the actual edges connection the cells (that the cells also reference)
        this.edgesInGrid.add(rightEdge);
      }
    } else {

      // ACCUMULATOR: Keeps track of all of the cells built for the maze
      ArrayList<ArrayList<ICell>> hiddenCells = new ArrayList<ArrayList<ICell>>();

      // ACCUMULATOR: Keeps track of all of the cells built for the maze
      // while keeping their identities. The cells in this list alias
      // those in _hiddenCells_ above in the same order, and are thereby
      // in a one-to-one correspondence
      ArrayList<ArrayList<Cell>> cells = new ArrayList<ArrayList<Cell>>();

      // Create and fill the grid of cells, initially disconnected
      // from each other but connected to 4 of their own `DummyCell`s
      for (int i = 0; i < cellsPerColumn; i += 1) {
        ArrayList<Cell> newRow = new ArrayList<Cell>();
        ArrayList<ICell> newHiddenRow = new ArrayList<ICell>();

        for (int n = 0; n < cellsPerRow; n += 1) {
          Cell newCell = new Cell();
          newRow.add(newCell);
          newHiddenRow.add(newCell);
        }
        cells.add(newRow);
        hiddenCells.add(newHiddenRow);
      }

      this.cells = hiddenCells;
      this.edgesInGrid = new ArrayList<CellEdge>();

      // Connect all cells not in the right or bottom edge to both their right and
      // bottom neighbors
      for (int i = 0; i < cellsPerColumn - 1; i += 1) {
        for (int j = 0; j < cellsPerRow - 1; j += 1) {

          Cell currCell = cells.get(i).get(j);
          Cell rightCell = cells.get(i).get(j + 1);
          Cell bottomCell = cells.get(i + 1).get(j);

          CellEdge rightEdge = currCell.connectAtRight(rightCell,
              this.weightMaker.nextHorizontalWeight());
          CellEdge bottomEdge = currCell.connectAtBottom(bottomCell,
              this.weightMaker.nextVerticalWeight());

          // Add these new edges to the edges in the maze. This edge aliases
          // the actual edge connecting the cells (that the cells also reference)
          this.edgesInGrid.add(rightEdge);
          this.edgesInGrid.add(bottomEdge);
        }
      }

      // Connect all cells in the bottom edge ONLY to their right neighbors
      for (int i = 0; i < cellsPerRow - 1; i += 1) {
        Cell currCell = cells.get(cellsPerColumn - 1).get(i);
        Cell rightCell = cells.get(cellsPerColumn - 1).get(i + 1);
        CellEdge rightEdge = currCell.connectAtRight(rightCell,
            this.weightMaker.nextHorizontalWeight());

        // Add the new edge to the edges in the maze. This edge aliases
        // the actual edges connection the cells (that the cells also reference)
        this.edgesInGrid.add(rightEdge);
      }

      // Connect all cells in the right edge ONLY to their bottom neighbors
      for (int i = 0; i < cellsPerColumn - 1; i += 1) {
        Cell currCell = cells.get(i).get(cellsPerRow - 1);
        Cell bottomCell = cells.get(i + 1).get(cellsPerRow - 1);
        CellEdge bottomEdge = currCell.connectAtBottom(bottomCell,
            this.weightMaker.nextVerticalWeight());
        this.edgesInGrid.add(bottomEdge);
      }
    }

    // Find all edges that will be open paths
    this.edgesInMST = new MazeUtils().allPathways(this.edgesInGrid,
        new ArrayUtils().foldl(this.cells, new ArrayList<ICell>(), new AppendList<ICell>()));

    this.setEndpoints();
  }

  // a convenience constructor for making mazes without consideration for bias but
  // considering hexagonal or square maze
  Maze(int cellsPerRow, int cellsPerColumn, int seed, boolean isHexagonalMaze) {
    this(cellsPerRow, cellsPerColumn, seed, false, false, isHexagonalMaze);
  }

  // a convenience constructor for making mazes without consideration for bias
  // (constructs
  // a square maze
  Maze(int cellsPerRow, int cellsPerColumn, int seed) {
    this(cellsPerRow, cellsPerColumn, seed, false, false, false);
  }

  // EFFECT: Sets one edge in the maze to be open and removes that edge from the
  // list of edges in the minimum search tree. If there are no more walls to
  // break down, the method does nothing
  void breakOneWall() {
    if (this.doneBreakingWalls()) {
      return;
    }

    this.edgesInMST.get(0).becomeEdgeInTree();
    this.edgesInMST.remove(0);
  }

  // EFFECT: Breaks all of the edges in the graph at once
  void breakAllWalls() {
    while (!this.doneBreakingWalls()) {
      this.breakOneWall();
    }
  }

  // Have we set every open pathway to the open state?
  public boolean doneBreakingWalls() {
    return this.edgesInMST.isEmpty();
  }

  // EFFECT: Set the starting and ending cells
  // to be displayed as their unique colors
  private void setEndpoints() {
    this.startingCell().enter(new StartState());
    this.finalCell().enter(new FinishState());
  }

  // NOTE: The following methods below (`drawInSceneHexagons()`,
  // `drawInSceneSquares()`, `drawInSceneInterpolatedHexagons()`,
  // and `drawInSceneInterpolatedSquares()` all have a similar form
  // but are NOT abstracted into a single helper method for performance
  // reasons. For grids above around 200x200 with 40_000 cells, the performance
  // cost of having another function call added up enough that animation speeds
  // slowed to an unacceptable level. We acknowledge that more efficient ways
  // of drawing the scene are likely to exist — our next strategy was to
  // incrementally draw
  // only those cells that needed drawing each frame and saving the previous scene
  // in `WorldMaze`. However, we did not have enough time to make the transition
  // with adequate testing in time for submission. This still works very well
  // for mazes on the order of 300x300.

  // EFFECT: Draws this Maze into the given scene of the given width and height
  // as if the cells were hexagons.This is much faster than creating
  // a large image representing the Maze with `BesideImage`s and `AboveImage`s.
  private void drawInSceneHexagons(WorldScene scene, int width, int height) {

    // Integer division intended (round down to fit the screen)
    // We choose the smaller of the two dimensions to ensure that the
    // cells fit in the screen. Since the height
    // of the hexagon is sqrt(3) * sideLength, we must also divide the WIDTH by
    // sqrt(3) value
    // to compensate (since the hexagon is rotated)
    final int normalizedWidth = (int) (Math.round(width / Math.sqrt(3)));
    final int normalizedHeight = 2 * height / 3; // Integer division intended

    // `effectiveCellsPerRow` is the number of cells in each row that
    // effectively take up space. Because the hexagons are offset
    // from each other, there is an additional half of a hexagon
    // of space that the last hexagon in the indented row juts out on the left
    // PLUS another half of a hexagon the un-indented row juts out on the right,
    // making
    // for one extra hexagon of space
    final int effectiveCellsPerRow = this.cellsPerRow + 1;

    // For the columns, because the each row is
    // in a honeycomb shape, the columns are actually compressed since
    // hexagon heights overlap. There are a total of
    // (n - 1) overlaps if there are n rows stacked on top of each other,
    // plus the additional half hexagon height jutting out at the top and
    // bottom to give (where h is the height of the hexagon and s is the side length
    // of the hexagon)
    //
    // h / 2 * (n - 1) + h / 2 + h / 2 =
    // h / 2 * (n - 1) + h =
    // s * (n - 1) + 2s = ns + s + 2s = ns + s
    // = (n + 1) s
    //
    // since h = 2s
    //
    // To see why, consider
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap + extra half jutting out
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap + extra half jutting out
    // the overlaps contribute half of the height of the hexagon, or simply `s`
    //
    // Thus, we have can fit
    //
    // s = (Height of scene) / (n + 1)
    final double effectiveCellsPerColumn = this.cellsPerColumn + 1;

    final int cellSideLength = (int) Math
        .round(Math.min(normalizedWidth / effectiveCellsPerRow, height / effectiveCellsPerColumn));

    // ACCUMULATOR: The X and Y locations in the scene to place the two given
    // images,
    // starting from the top left corner and working our way down the rows
    // of the cells
    double x = 0;
    double y = 0;

    // The amount by which we should move each time we add a new hexagon
    // in both the X and Y directions
    final double deltaX = Math.round(Math.sqrt(3) * cellSideLength) - 1;
    final double deltaY = 3 * cellSideLength / 2 - 1;

    // ACCUMULATOR: Keeps track of whether or not we are located
    // in a row that is offset with respect to neighboring rows
    boolean inOffsetRow = false;

    for (ArrayList<ICell> cellsInRow : this.cells) {
      for (ICell cell : cellsInRow) {

        // Move the pinhole such that the top-left corner
        // of the cell is the pinhole
        WorldImage drawnCell = cell.toImage(cellSideLength);
        WorldImage cellImage = drawnCell.movePinhole(-drawnCell.getWidth() / 2,
            -drawnCell.getHeight() / 2);

        scene.placeImageXY(cellImage, (int) (Math.round(x)), (int) (Math.round(y)));

        // Move over the sidelength of the cell
        x += deltaX;
      }

      if (inOffsetRow) {
        // Move back to the beginning of the row
        x = 0;
      } else {
        // In this case, we are moving to a row
        // that is offset from this row. Hence, we start at
        // HALF deltaX (that is, a single height of one of the equilateral
        // triangles in the hexagon)
        x = deltaX / 2;
      }

      // Move down to the next row
      y += deltaY;

      // We are now in or not in an offset row
      inOffsetRow = !inOffsetRow;
    }
  }

  // EFFECT: Draws this Maze into the given scene of the given width and height
  // as if the cells were squares.This is much faster than creating
  // a large image representing the Maze with `BesideImage`s and `AboveImage`s.
  private void drawInSceneSquares(WorldScene scene, int width, int height) {
    // Integer division intended (round down to fit the screen)
    // We choose the smaller of the two dimensions to ensure that the
    // cells fit in the screen
    final int cellSideLength = Math.min(width / this.cellsPerRow, height / this.cellsPerColumn);

    // ACCUMULATOR: The X and Y locations in the scene to place the two given
    // images,
    // starting from the top left corner and working our way down the rows
    // of the cells
    int x = 0;
    int y = 0;

    // Draw each row of cells and stack them on top of each other
    for (ArrayList<ICell> cellsInRow : this.cells) {
      for (ICell cell : cellsInRow) {

        // Move the pinhole such that the top-left corner
        // of the cell is the pinhole
        WorldImage drawnCell = cell.toImage(cellSideLength);
        WorldImage cellImage = drawnCell.movePinhole(-drawnCell.getWidth() / 2,
            -drawnCell.getHeight() / 2);

        scene.placeImageXY(cellImage, x, y);

        // Move over the sidelength of the cell
        x += cellSideLength;
      }

      // Move down the side length of the cell and move
      // back to the beginning of the row
      x = 0;
      y += cellSideLength;
    }
  }

  // EFFECT: Draws this Maze into the given scene of the given width and height.
  // This is much faster than creating
  // a large image representing the Maze with `BesideImage`s and `AboveImage`s.
  void drawInScene(WorldScene scene, int width, int height) {
    if (this.isHexagonalMaze) {
      this.drawInSceneHexagons(scene, width, height);
    } else {
      this.drawInSceneSquares(scene, width, height);
    }
  }

  // EFFECT: Draws this Maze into the given scene using interpolated
  // values for cells by computing their distances from the start or end. This is
  // much faster than
  // creating a large image representing the Maze with `BesideImage`s and
  // `AboveImage`s
  void drawInSceneInterpolated(WorldScene scene, int width, int height, boolean fromStart) {
    if (this.isHexagonalMaze) {
      this.drawInSceneInterpolatedHexagons(scene, width, height, fromStart);
    } else {
      this.drawInSceneInterpolatedSquares(scene, width, height, fromStart);
    }
  }

  // EFFECT: Draws this Maze into the given scene using interpolated
  // values for hexagonal cells by computing their distances from the start or
  // end. This is much faster than
  // creating a large image representing the Maze with `BesideImage`s and
  // `AboveImage`s
  private void drawInSceneInterpolatedHexagons(WorldScene scene, int width, int height,
      boolean fromStart) {
    // we declare the snapshot before assignment to ensure it's in the scope of
    // the remainder of the method. It will be assigned, as there is an if else that
    // assigns it to a snapshot in both cases.
    ExhaustedSnapshot snapshotThroughEntireMaze;

    if (fromStart) {
      snapshotThroughEntireMaze = this.newExhaustedSnapshotFromStart();
    } else {
      snapshotThroughEntireMaze = this.newExhaustedSnapshotFromFinish();
    }

    // Integer division intended (round down to fit the screen)
    // We choose the smaller of the two dimensions to ensure that the
    // cells fit in the screen. Since the height
    // of the hexagon is sqrt(3) * sideLength, we must also divide the WIDTH by
    // sqrt(3) value
    // to compensate (since the hexagon is rotated)
    final int normalizedWidth = (int) (Math.round(width / Math.sqrt(3)));
    final int normalizedHeight = 2 * height / 3; // Integer division intended

    // `effectiveCellsPerRow` is the number of cells in each row that
    // effectively take up space. Because the hexagons are offset
    // from each other, there is an additional half of a hexagon
    // of space that the last hexagon in the indented row juts out on the left
    // PLUS another half of a hexagon the un-indented row juts out on the right,
    // making
    // for one extra hexagon of space
    final int effectiveCellsPerRow = this.cellsPerRow + 1;

    // For the columns, because the each row is
    // in a honeycomb shape, the columns are actually compressed since
    // hexagon heights overlap. There are a total of
    // (n - 1) overlaps if there are n rows stacked on top of each other,
    // plus the additional half hexagon height jutting out at the top and
    // bottom to give (where h is the height of the hexagon and s is the side length
    // of the hexagon)
    //
    // h / 2 * (n - 1) + h / 2 + h / 2 =
    // h / 2 * (n - 1) + h =
    // s * (n - 1) + 2s = ns + s + 2s = ns + s
    // = (n + 1) s
    //
    // since h = 2s
    //
    // To see why, consider
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap + extra half jutting out
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap
    // ⬡⬡⬡⬡⬡⬡⬡⬡ -> overlap + extra half jutting out
    // the overlaps contribute half of the height of the hexagon, or simply `s`
    //
    // Thus, we have can fit
    //
    // s = (Height of scene) / (n + 1)
    final double effectiveCellsPerColumn = this.cellsPerColumn + 1;

    final int cellSideLength = (int) Math
        .round(Math.min(normalizedWidth / effectiveCellsPerRow, height / effectiveCellsPerColumn));

    // ACCUMULATOR: The X and Y locations in the scene to place the two given
    // images,
    // starting from the top left corner and working our way down the rows
    // of the cells
    double x = 0;
    double y = 0;

    // The amount by which we should move each time we add a new hexagon
    // in both the X and Y directions
    final double deltaX = Math.round(Math.sqrt(3) * cellSideLength) - 1.5;
    final double deltaY = 3 * cellSideLength / 2 - 1;

    // ACCUMULATOR: Keeps track of whether or not we are located
    // in a row that is offset with respect to neighboring rows
    boolean inOffsetRow = false;

    for (ArrayList<ICell> cellsInRow : this.cells) {
      for (ICell cell : cellsInRow) {

        // Move the pinhole such that the top-left corner
        // of the cell is the pinhole
        WorldImage drawnCell = cell.toImage(cellSideLength,
            snapshotThroughEntireMaze.distanceFromStart(cell),
            snapshotThroughEntireMaze.distanceToFarthest());
        WorldImage cellImage = drawnCell.movePinhole(-drawnCell.getWidth() / 2,
            -drawnCell.getHeight() / 2);

        scene.placeImageXY(cellImage, (int) (Math.round(x)), (int) (Math.round(y)));

        // Move over the side length of the cell
        x += deltaX;
      }

      if (inOffsetRow) {
        // Move back to the beginning of the row
        x = 0;
      } else {
        // In this case, we are moving to a row
        // that is offset from this row. Hence, we start at
        // HALF deltaX (that is, a single height of one of the equilateral
        // triangles in the hexagon)
        x = deltaX / 2;
      }

      // Move down to the next row
      y += deltaY;

      // We are now in or not in an offset row
      inOffsetRow = !inOffsetRow;
    }
  }

  // EFFECT: Draws this Maze into the given scene of the given width and height
  // as if the cells were squares.This is much faster than creating
  // a large image representing the Maze with `BesideImage`s and `AboveImage`s.
  private void drawInSceneInterpolatedSquares(WorldScene scene, int width, int height,
      boolean fromStart) {

    // we declare the snapshot before assignment to ensure it's in the scope of
    // the remainder of the method. It will be assigned, as there is an if else that
    // assigns it to a snapshot in both cases.
    ExhaustedSnapshot snapshotThroughEntireMaze;

    if (fromStart) {
      snapshotThroughEntireMaze = this.newExhaustedSnapshotFromStart();
    } else {
      snapshotThroughEntireMaze = this.newExhaustedSnapshotFromFinish();
    }

    // Integer division intended (round down to fit the screen)
    // We choose the smaller of the two dimensions to ensure that the
    // cells fit in the screen
    final int cellSideLength = Math.min(width / this.cellsPerRow, height / this.cellsPerColumn);

    // ACCUMULATOR: The X and Y locations in the scene to place the two given
    // images,
    // starting from the top left corner and working our way down the rows
    // of the cells
    int x = 0;
    int y = 0;

    // Draw each row of cells and stack them on top of each other
    for (ArrayList<ICell> cellsInRow : this.cells) {
      for (ICell cell : cellsInRow) {

        // Move the pinhole such that the top-left corner
        // of the cell is the pinhole
        WorldImage drawnCell = cell.toImage(cellSideLength,
            snapshotThroughEntireMaze.distanceFromStart(cell),
            snapshotThroughEntireMaze.distanceToFarthest());
        WorldImage cellImage = drawnCell.movePinhole(-drawnCell.getWidth() / 2,
            -drawnCell.getHeight() / 2);

        scene.placeImageXY(cellImage, x, y);

        // Move over the side length of the cell
        x += cellSideLength;
      }

      // Move down the side length of the cell and move
      // back to the beginning of the row
      x = 0;
      y += cellSideLength;
    }
  }

  // Produces the starting DepthFirstSnapshot for this maze
  DepthFirstSnapshot newDepthFirstSnapshot() {
    return new DepthFirstSnapshot(this.startingCell(), this.finalCell());
  }

  // Produces the starting BreadthFirstSnapshot for this maze
  BreadthFirstSnapshot newBreadthFirstSnapshot() {
    return new BreadthFirstSnapshot(this.startingCell(), this.finalCell());
  }

  // Produces an exhaustive snapshot through this entire maze
  // starting from the first cell in the maze
  ExhaustedSnapshot newExhaustedSnapshotFromStart() {
    return new ExhaustedSnapshot(this.startingCell());
  }

  // Produces an exhaustive snapshot through this entire maze
  // starting from the last cell in the maze
  ExhaustedSnapshot newExhaustedSnapshotFromFinish() {
    return new ExhaustedSnapshot(this.finalCell());
  }

  // Produces the starting ManualSnapshot for this maze
  ManualSearchSnapshot newManualSnapshot() {
    return new ManualSearchSnapshot(this.startingCell(), this.finalCell());
  }

  // Produces the starting cell in the maze
  private ICell startingCell() {
    return this.cells.get(0).get(0);
  }

  // Produces the goal/final cell in the maze
  private ICell finalCell() {
    return this.cells.get(this.cellsPerColumn - 1).get(this.cellsPerRow - 1);
  }
}

// Represents a utility class for maze
final class MazeUtils {

  // Given all of the individual cells in a maze, as well as all of the
  // connections
  // between individual cells in the maze, determines that subset of edges which
  // form a minimum spanning tree through the cells, encoded as a list without any
  // particular ordering to the edges
  ArrayList<CellEdge> allPathways(ArrayList<CellEdge> allEdges, ArrayList<ICell> allCells) {

    ArrayList<CellEdge> edgesInMST = new ArrayList<CellEdge>();

    // Make a copy of the edges in the Maze connecting the cells and sort it
    ArrayList<CellEdge> sortedEdges = new ArrayList<CellEdge>(allEdges);
    sortedEdges.sort(new ByEdgeWeight());

    // Now create a UnionFind data structure which uses the cells as nodes
    UnionFind<ICell> cellForest = new UnionFind<ICell>(allCells);

    // Now, process the edges one-by-one
    for (CellEdge edge : sortedEdges) {

      ICell rootFirst = cellForest.findRepresentative(edge.first);
      ICell rootSecond = cellForest.findRepresentative(edge.second);

      // If the cells in this edge are not already in the same
      // trees in the UnionFind forest, we know adding this edge
      // won't create a cycle
      if (rootFirst != rootSecond) {
        edgesInMST.add(edge);
        cellForest.union(rootFirst, rootSecond);
      }
    }
    return edgesInMST;
  }
}

// A `UnionFind<T>` is a representation of a forest of trees
// each of whose nodes are mutually disjoint from one another that
// supports efficient access for locating a node's representative
// for its respective tree and for combining trees together in the forest.
final class UnionFind<T> {

  // This is an encoding of the forest represented by
  // this `UnionFind`. Each node in the forest is identified
  // by its intentional identity
  private final IdentityHashMap<T, T> nodeLinks;

  // This is a mapping of distinct objects to the
  // number of other items in the forest encoded by the
  // `nodeLinks` structure that refer to the given node.
  // Note, though, that values mapped in this structure to given instances are not
  // necessarily meaningful UNLESS they are current representatives
  // for the forests in the tree (see the discussion for
  private final IdentityHashMap<T, Integer> treeSizes;

  // Construct a new `UnionFind` which builds trees from the elements
  // of the given list. The starting forest is a collection of single-node
  // trees, one for each item in the list of elements
  UnionFind(ArrayList<T> startList) {

    IdentityHashMap<T, T> nodeLinks = new IdentityHashMap<T, T>();
    IdentityHashMap<T, Integer> treeSizes = new IdentityHashMap<T, Integer>();

    for (T element : startList) {
      nodeLinks.put(element, element);
      treeSizes.put(element, 1);
    }

    this.nodeLinks = nodeLinks;
    this.treeSizes = treeSizes;
  }

  // Is the given item contained in the set of nodes visible to this UnionFind?
  boolean contains(T element) {
    return this.nodeLinks.containsKey(element);
  }

  // Returns the top of the tree that the given node is a part of
  T findRepresentative(T node) {

    // ensure that no item can be referenced if it did not begin in the original
    // node links.
    // Nodes without representatives can only occur at the first reference to that
    // item,
    // and not farther up any trees because we will never use an item that is not in
    // the
    // `nodeLinks` to go beyond the original node. That is, since we ensure
    if (!this.contains(node)) {
      throw new IllegalArgumentException("Could not locate a representative for the"
          + " provided item as it is not in the set of items visible to this UnionFind instance");
    }

    T current = node;
    while (this.nodeLinks.get(current) != current) {
      current = this.nodeLinks.get(current);
    }
    return current;
  }

  // EFFECT: Given two elements in the forest of items
  // described by this UnionFind, combines the two trees
  // the two elements are apart of by merging the smaller
  // of the two trees into the larger of the two
  void union(T element1, T element2) {

    // Again, ensure that both items are within the set
    // of elements visible to this UnionFind
    if (!this.contains(element1)) {
      throw new IllegalArgumentException("Could not perform a union for the"
          + " first provided item as it is not in the set of items "
          + "visible to this UnionFind instance");
    }

    if (!this.contains(element2)) {
      throw new IllegalArgumentException("Could not perform a union for the"
          + " second provided item as it is not in the set of items "
          + "visible to this UnionFind instance");
    }

    // First, determine the representatives of the two elements in the forest
    T root1 = this.findRepresentative(element1);
    T root2 = this.findRepresentative(element2);

    // Then, determine which representative has a larger
    // cohort (number of nodes that refer to it as its representative)
    int numNodesInRoot1Tree = this.treeSizes.get(root1);
    int numNodesInRoot2Tree = this.treeSizes.get(root2);

    // Absorb the smaller of the two trees into the larger of the two
    // by having the appropriate root refer to the other root as its
    // new representative. Also, add the nodes in the smaller
    // of the two trees to the count of the larger of the two trees.
    // If the trees have the same size, the second element becomes a child
    // of the first
    if (numNodesInRoot1Tree < numNodesInRoot2Tree) {
      this.nodeLinks.put(root1, root2);
      this.treeSizes.put(root2, numNodesInRoot1Tree + numNodesInRoot2Tree);
    } else {
      this.nodeLinks.put(root2, root1);
      this.treeSizes.put(root1, numNodesInRoot1Tree + numNodesInRoot2Tree);
    }
  }
}

// An `AddFind<T>` is a data structure that represents
// a rooted tree with edges pointing away from the root, also known
// as an anti-arborescence or in-tree (see https://en.wikipedia.org/wiki/Tree_(graph_theory), subsection
// "Rooted Tree" for more details). An `AddFind<T>` efficiently supports adding new
// items to the in-tree as well as for finding the depth of a given node in the tree, as
// well as finding the ancestors of a node in the in-tree. The structure is related
// to the UnionFind except that is designed to be an extensible tree rather than representing
// a forest from a fixed set of nodes
final class AddFind<T> implements Iterable<T> {

  // An encoding of the tree structures that this AddFind represents.
  // Items are mapped to their parents using their intentional identity
  // and the root is mapped to itself
  private final IdentityHashMap<T, T> nodePaths;

  // An encoding of tree heights of each node in the tree.
  // Each node is mapped to its depth in the in-tree.
  private final IdentityHashMap<T, Integer> nodeDepths;

  // Create a new AddFind without any nodes in the tree
  AddFind() {
    this.nodePaths = new IdentityHashMap<T, T>();
    this.nodeDepths = new IdentityHashMap<T, Integer>();
  }

  // EFFECT: Adds an element to the AddFind, with _root_ being added as
  // the first element in the AddFind
  //
  // Note: This can only be used to add the first element into an AddFind,
  // and should be used exclusively to add the first element into an AddFind,
  // as the given item will be the root to all other items
  void addRoot(T root) {
    if (this.isEmpty()) {
      this.nodePaths.put(root, root);

      // The depth of the root is 0
      this.nodeDepths.put(root, 0);
    } else {
      throw new UnsupportedOperationException("Cannot add a root to a non-empty AddFind");
    }
  }

  // EFFECT: Adds an element to the AddFind, with _element_ being added
  // as a child of _parent_ if _parent_ is already in the AddFind. Otherwise,
  // an exception is thrown if _parent_ is NOT part of the add find or if
  // the element is intentionally equal to the element we are adding
  void add(T element, T parent) {
    if (element == parent) {
      throw new IllegalArgumentException("Cannot directly add an item to the `AddFind<T>` "
          + "whose parent is itself since that element is not the root");
    }

    // Ensure duplicate items are not added to the AddFind
    if (this.contains(element)) {
      throw new IllegalArgumentException(
          "Cannot add the given element to the AddFind because it " + "is already in the AddFind");
    }

    if (this.contains(parent)) {
      this.nodePaths.put(element, parent);

      // The depth of this new element is one more than its parent
      this.nodeDepths.put(element, this.nodeDepths.get(parent) + 1);
    } else {
      throw new IllegalArgumentException("Cannot add an element to the AddFind"
          + " to a node that is not already part of the AddFind");
    }
  }

  // Returns how far down in the AddFind in-tree that the given node resides
  int depthOf(T element) {
    if (!this.contains(element)) {
      throw new IllegalArgumentException(
          "The given element is not in the AddFind to get a depth for");
    }

    return this.nodeDepths.get(element);
  }

  // Are there any elements in this AddFind?
  boolean isEmpty() {
    return this.nodePaths.isEmpty();
  }

  // Is the given item contained in this AddFind?
  boolean contains(T element) {
    return this.nodePaths.containsKey(element);
  }

  // Returns an ArrayList of the complete path of parents from a
  // given item to the starting item
  ArrayList<T> getLineage(T element) {

    if (!this.contains(element)) {
      throw new IllegalArgumentException("The given element is not present in the AddFind");
    }

    // ACCUMULATOR: The path of ancestors
    // that the given element has to the root
    ArrayList<T> pathBack = new ArrayList<T>();

    T currentElement = element;

    while (currentElement != this.nodePaths.get(currentElement)) {
      pathBack.add(currentElement);
      currentElement = this.nodePaths.get(currentElement);
    }

    // We use the root as our stopping point, so the root
    // is not added in the `while` loop. So we add it here
    // as it is one of the ancestors of _element_
    pathBack.add(currentElement);

    return pathBack;
  }

  // Creates an iterator over the elements in this AddFind
  // Note: The order of traversal of the iterator is not
  // guaranteed, and should only be used to get access to all
  // nodes in the AddFind regardless of order
  public Iterator<T> iterator() {
    return this.nodePaths.keySet().iterator();
  }
}

// Represents the status of a node search algorithm as it
// would progress
abstract class ASearchSnapshot {

  // All of the paths that we have looked through under this algorithm
  protected AddFind<ICell> travelledPaths;

  // The cells we still need to process
  protected ICollection<ICell> workList;

  private int maxDistanceFromStart;

  protected ICell start;
  protected ICell finish;
  protected ICell current;

  ASearchSnapshot(ICell start, ICell finish, ICollection<ICell> collection) {
    this.travelledPaths = new AddFind<ICell>();
    this.workList = collection;
    this.workList.add(start);
    this.maxDistanceFromStart = 0;
    this.start = start;
    this.finish = finish;
    this.current = start;
  }

  // Returns whether or not the search has reached the designated final cell
  // EFFECT: Pulls an item off of the _workList_, adds the item's viable
  // neighbors,
  // and records seeing this item
  boolean searchOneStep() {

    // find the item closest to removal from the collection that is a part of
    // an unseen path
    this.current = this.nextItem();

    // If, for some reason, _current_ has already been seen before in the
    // AddFind (that is, if the algorithm has processed this node), then
    // we know that the promise kept by the users of this class was broken, and we
    // have
    // inadvertently hit a cycle in the ICell graph when none were expected
    if (this.travelledPaths.contains(current)) {
      throw new IllegalStateException("The search snapshot's algorithm reached a cycle in "
          + "the graph it is traversing when none were expected");
    }

    if (this.travelledPaths.isEmpty()) {
      this.travelledPaths.addRoot(this.current);
    } else {
      // find all open paths extending from this cell
      ArrayList<ICell> outCells = this.current.adjacentOutCells();

      // determine the cell that placed this cell into the collection by finding
      // which of its neighbors have been seen
      //
      // In a minimum spanning tree, there is only one path connecting two items,
      // so we can be sure that only one of _current_'s neighbors have been seen.
      // Since only one neighbor to _current_ has been seen, that has to be the
      // cell that placed _current_ into the collection.
      ICell parent = new ArrayUtils().firstToSatisfy(outCells,
          new AddFindContains<ICell>(this.travelledPaths));

      // Add this pathway to our AddFind
      this.travelledPaths.add(this.current, parent);
    }

    this.maxDistanceFromStart = Math.max(this.maxDistanceFromStart,
        this.travelledPaths.depthOf(this.current));

    // Taking all potential pathways from _current_, we will remove the parent
    // cell of this cell and push the remaining neighboring cells into the
    // collection
    ArrayList<ICell> nextToAdd = new ArrayUtils().filter(this.current.adjacentOutCells(),
        new NotPred<ICell>(new AddFindContains<ICell>(this.travelledPaths)));

    for (ICell cell : nextToAdd) {
      this.workList.add(cell);
    }

    // Have we finished the search (have we found the finish?)?
    return this.current == this.finish;
  }

  // Are there still more items left to search for? (is our _workList_ empty?)
  boolean hasMoreToSee() {
    return !this.workList.isEmpty();
  }

  // determines how far the farthest point from the start is
  final int distanceToFarthest() {
    return this.maxDistanceFromStart;
  }

  // Returns the next item on the stack that has not yet been seen
  // EFFECT: removes items off the stack until finding a new item
  private ICell nextItem() {

    if (!this.hasMoreToSee()) {
      throw new IllegalStateException(
          "The worklist is empty, more items cannot" + " be popped off.");
    }

    // take an item off the top of the stack
    ICell maybeCurrent = this.workList.remove();

    // if we've already seen this item, pop a new item off the stack
    // until finding an item that is a part of a new path
    while (this.travelledPaths.contains(maybeCurrent)) {

      if (!this.hasMoreToSee()) {
        throw new IllegalStateException(
            "The worklist is empty, more items cannot" + " be popped off.");
      }
      maybeCurrent = this.workList.remove();
    }

    return maybeCurrent;
  }

  // EFFECT: Marks all cells that have been seen with the Seen state, and
  // marks all heads of new paths with the Path state
  void makeAllSeen() {

    // produce an iterator over all items seen from our AddFind
    Iterator<ICell> seenIter = this.travelledPaths.iterator();

    // for all cells seen, set state to seen
    while (seenIter.hasNext()) {
      ICell curr = seenIter.next();
      curr.enter(new SeenState());
    }

    // for our current cell, set state to head
    this.current.enter(new HeadState());
  }

  // EFFECT: marks just the head as a head state
  void makeOnlyHead() {

    // produce an iterator over all items seen from our AddFind
    Iterator<ICell> seenIter = this.travelledPaths.iterator();

    // for all cells seen, set state to unseen
    while (seenIter.hasNext()) {
      ICell curr = seenIter.next();
      curr.enter(new UnseenState());
    }

    // for our current cell, set state to head
    this.current.enter(new HeadState());
  }

  // EFFECT: Marks all cells that have been seen with the Seen state, and
  // marks all heads of new paths with the Path state
  void makeAllUnseen() {
    // Produce an iterator over all items seen from our AddFind
    Iterator<ICell> seenIter = this.travelledPaths.iterator();

    // for all cells seen, query the cells to let them
    // know they're no longer seen
    while (seenIter.hasNext()) {
      ICell curr = seenIter.next();
      curr.enter(new UnseenState());
    }

    // Reset the starting and ending cells
    this.start.enter(new StartState());
    this.finish.enter(new FinishState());
  }

  // EFFECT: Finds all items that lead to the final item and gives each item
  // in that list the Seen State
  void makeAllInPath() {
    ArrayList<ICell> finishedPath = this.travelledPaths.getLineage(finish);

    for (ICell cell : finishedPath) {
      cell.enter(new PathState());
    }
  }

  // Returns the distance of the given cell from this ASearchSnapshot's starting
  // cell.
  // If the cell has not been visited by the search algorithm described by this
  // object,
  // an exception is thrown
  final int distanceFromStart(ICell cell) {
    if (!this.travelledPaths.contains(cell)) {
      throw new IllegalArgumentException(
          "The search algorithm state has not yet processed the given ICell");
    }

    return this.travelledPaths.depthOf(cell);
  }

  // determines if we can move up, and will move up or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveUp() {
    throw new UnsupportedOperationException("Cannot move up during a automatic search.");
  }

  // determines if we can move down, and will move down or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveDown() {
    throw new UnsupportedOperationException("Cannot move down during a automatic search.");
  }

  // determines if we can move left, and will move left or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveLeft() {
    throw new UnsupportedOperationException("Cannot move left during a automatic search.");
  }

  // determines if we can move right, and will move right or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveRight() {
    throw new UnsupportedOperationException("Cannot move right during a automatic search.");
  }

  // determines if we can move top left, and will move top left or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveTopLeft() {
    throw new UnsupportedOperationException("Cannot move top left during a automatic search.");
  }

  // determines if we can move top right, and will move top right or not respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveTopRight() {
    throw new UnsupportedOperationException("Cannot move top right during a automatic search.");
  }

  // determines if we can move bottom left, and will move bottom left or not
  // respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveBottomLeft() {
    throw new UnsupportedOperationException("Cannot move bottom left during a automatic search.");
  }

  // determines if we can move bottom right, and will move bottom right or not
  // respond
  // NOTE: for all searches except for manual, we do not want this to
  // be called, and will throw an exception accordingly
  public void maybeMoveBottomRight() {
    throw new UnsupportedOperationException("Cannot move bottomRight during a automatic search.");
  }
}

//Represents the status of a manual search as it would progress
class ManualSearchSnapshot extends ASearchSnapshot {

  ManualSearchSnapshot(ICell start, ICell finish) {
    // we choose to pass in a stack arbitrarily, we will not use the
    // stack but require a collection for our abstract class, which
    // we still intend to use
    super(start, finish, new OurStack<ICell>());
  }

  // Are there still more items left to search for?
  // In manual search, you will not run into an issue of an empty worklist
  // preventing moving, so we will always be able to see more
  public boolean hasMoreToSee() {
    return true;
  }

  // Returns whether or not the search has reached the designated final cell
  // EFFECT: No effect for manual search
  public boolean searchOneStep() {

    return this.current == this.finish;
  }

  // moves to the given cell if it is different from the current cell, not moving
  // if the cells are the same
  public void maybeMoveTo(ICell next) {

    // if these are the same, up is an illegal path, and we do not move
    if (next == this.current) {
      return;
    }

    if (!this.travelledPaths.contains(next)) {

      if (this.travelledPaths.isEmpty()) {
        this.travelledPaths.addRoot(this.current);
      }

      this.travelledPaths.add(next, this.current);

    }

    this.current = next;
  }

  // EFFECT: determines if we can move up, and will move up or not respond
  public void maybeMoveUp() {

    this.maybeMoveTo(this.current.maybeTop());
  }

  // determines if we can move down, and will move down or not respond
  public void maybeMoveDown() {
    this.maybeMoveTo(this.current.maybeBottom());
  }

  // determines if we can move left, and will move left or not respond
  public void maybeMoveLeft() {
    this.maybeMoveTo(this.current.maybeLeft());
  }

  // determines if we can move right, and will move right or not respond
  public void maybeMoveRight() {
    this.maybeMoveTo(this.current.maybeRight());
  }

  // determines if we can move top left, and will move top left or not respond
  public void maybeMoveTopLeft() {
    this.maybeMoveTo(this.current.maybeTopLeft());
  }

  // determines if we can move top right, and will move top right or not respond
  public void maybeMoveTopRight() {
    this.maybeMoveTo(this.current.maybeTopRight());
  }

  // determines if we can move bottom left, and will move bottom left or not
  // respond
  public void maybeMoveBottomLeft() {
    this.maybeMoveTo(this.current.maybeBottomLeft());
  }

  // determines if we can move bottom right, and will move bottom right or not
  // respond
  public void maybeMoveBottomRight() {
    this.maybeMoveTo(this.current.maybeBottomRight());
  }
}

// Represents the status of a depth-first node search algorithm as it would progress
class DepthFirstSnapshot extends ASearchSnapshot {
  DepthFirstSnapshot(ICell start, ICell finish) {
    super(start, finish, new OurStack<ICell>());
  }
}

// Represents a completely exhausted depth-first node search algorithm
// snapshot as it progressed through the entire maze
class ExhaustedSnapshot extends DepthFirstSnapshot {
  ExhaustedSnapshot(ICell start) {
    super(start, new Cell());

    // Keep going until we've seen everything in the structure
    // or until we've hit the final node that we care about.
    // EFFECT: this.searchOneStep has the side effect of
    // progressing the search forward
    // TERMINATION: The cells are part of a minimum spanning tree
    while (this.hasMoreToSee() && !this.searchOneStep()) {
    }

  }
}

// Represents the status of a breadth-first node search algorithm as it would progress
class BreadthFirstSnapshot extends ASearchSnapshot {

  BreadthFirstSnapshot(ICell start, ICell finish) {
    super(start, finish, new Queue<ICell>());
  }

  // EFFECT: Marks all cells that have been seen with the Seen state, and
  // marks all heads of new paths with the Path state
  void makeAllSeen() {
    super.makeAllSeen();

    for (ICell cell : this.workList) {
      cell.enter(new FutureHeadState());
    }
  }

  // EFFECT: marks just the head as a head state
  void makeOnlyHead() {

    super.makeOnlyHead();

    for (ICell cell : this.workList) {
      cell.enter(new UnseenState());
    }
  }

  // EFFECT: Marks all cells that have been seen with the Seen state, and
  // marks all heads of new paths with the Path state
  void makeAllUnseen() {
    super.makeAllUnseen();
    for (ICell cell : this.workList) {
      cell.enter(new UnseenState());
    }
  }
}

// Represents a mutable collection of items
interface ICollection<T> extends Iterable<T> {
  // Is this collection empty?
  boolean isEmpty();

  // EFFECT: Adds the item to the collection
  void add(T item);

  // Returns the first item of the collection
  // EFFECT: Removes that first item
  T remove();
}

// Represents a version of a stack that implements ICollection
class OurStack<T> implements ICollection<T> {

  private Stack<T> contents;

  OurStack() {
    this.contents = new Stack<T>();
  }

  // Is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  // Returns the first item of the collection
  // EFFECT: removes that first item
  public T remove() {

    if (this.isEmpty()) {
      throw new EmptyStackException();
    }

    return this.contents.pop();
  }

  // EFFECT: adds the item to the collection
  public void add(T item) {
    this.contents.push(item);
  }

  // produces a iterator over our stack
  public Iterator<T> iterator() {
    return this.contents.iterator();
  }
}

// Represents a version of a queue that implements ICollection
class Queue<T> implements ICollection<T> {
  private Deque<T> contents;

  Queue() {
    this.contents = new ArrayDeque<T>();
  }

  // Is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  // Returns the first item of the collection
  // EFFECT: removes that first item
  public T remove() {
    if (this.isEmpty()) {
      throw new NoSuchElementException();
    }
    return this.contents.removeFirst();
  }

  // EFFECT: adds the item to the collection
  public void add(T item) {
    this.contents.addLast(item);
  }

  // produces an iterator over our Queue
  public Iterator<T> iterator() {
    return this.contents.iterator();
  }
}

// Represents a single unit in a maze
interface ICell {
  // returns a list of cells that are connected to this cell by open edges
  ArrayList<ICell> adjacentOutCells();

  // returns a list of cells that are connected to this cell
  ArrayList<ICell> allAdjacentCells();

  // EFFECT: Changes the state of this cell, reassigning it to the given state
  void enter(ICellState state);

  // Renders this cell as image that fits in a square of the given side length
  WorldImage toImage(int sideLength);

  // Renders this cell as image that fits in a square of the given side length
  // and as if it were the given distance away from some starting point with
  // respect to some maximum distance from the given point
  WorldImage toImage(int sideLength, int distanceFromLoc, int maxDistanceFromLoc);

  // returns the cell above this cell, if there is an open connection to it,
  // returning this cell otherwise
  ICell maybeTop();

  // returns the cell below this cell, if there is an open connection to it,
  // returning this cell otherwise
  ICell maybeBottom();

  // returns the cell to the left of this cell, if there is an open connection to
  // it,
  // returning this cell otherwise
  ICell maybeLeft();

  // returns the cell to the right of this cell, if there is an open connection to
  // it,
  // returning this cell otherwise
  ICell maybeRight();

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  ICell maybeTopRight();

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  ICell maybeBottomRight();

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  ICell maybeBottomLeft();

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  ICell maybeTopLeft();
}

//represents a square ICell, with junctions on each of the four edges,
//connecting it to it's neighbors. Cells on the maze's boundary are connected
//to dummy cells that form a grid boundary.
class Cell implements ICell {
  private CellEdge top;
  private CellEdge left;
  private CellEdge right;
  private CellEdge bottom;
  private ICellState state;

  // constructs a cell with connections to four dummy cells
  Cell() {
    this.top = new CellEdge(this, new DummyCell());
    this.left = new CellEdge(this, new DummyCell());
    this.right = new CellEdge(this, new DummyCell());
    this.bottom = new CellEdge(this, new DummyCell());
    this.state = new UnseenState();
  }

  // Returns the new right connection between this cell and the given one
  // EFFECT: Connects this cell with the given one on its right
  CellEdge connectAtRight(Cell rightCell, int weight) {
    CellEdge connection = new CellEdge(this, rightCell, weight);
    this.right = connection;
    rightCell.left = connection;
    return connection;
  }

  // Returns the new right connection between this cell and the given one
  // EFFECT: Connects this cell with the given one on its right
  CellEdge connectAtBottom(Cell bottomCell, int weight) {
    CellEdge connection = new CellEdge(this, bottomCell, weight);
    this.bottom = connection;
    bottomCell.top = connection;
    return connection;
  }

  // Renders this Cell with the given wall color and background color
  private WorldImage toImage(int sideLength, Color cellColor, Color wallColor) {
    WorldImage drawnCell = new RectangleImage(sideLength, sideLength, OutlineMode.SOLID, cellColor);

    if (!this.top.inTree()) {
      WorldImage wall = new LineImage(new Posn(sideLength, 0), wallColor);
      drawnCell = new OverlayOffsetImage(wall, 0, (double) (sideLength) / 2.0, drawnCell);
    }
    if (!this.left.inTree()) {
      WorldImage wall = new LineImage(new Posn(0, sideLength), wallColor);
      drawnCell = new OverlayOffsetImage(wall, (double) (sideLength) / 2.0, 0, drawnCell);
    }
    if (!this.right.inTree()) {
      WorldImage wall = new LineImage(new Posn(0, sideLength), wallColor);
      drawnCell = new OverlayOffsetImage(wall, (double) (-sideLength) / 2.0, 0, drawnCell);
    }
    if (!this.bottom.inTree()) {
      WorldImage wall = new LineImage(new Posn(sideLength, 0), wallColor);
      drawnCell = new OverlayOffsetImage(wall, 0, (double) (-sideLength) / 2.0, drawnCell);
    }
    return drawnCell;
  }

  // renders this cell as an image
  public WorldImage toImage(int sideLength) {
    return this.toImage(sideLength, this.state.standardColor(), Constants.WALL_COLOR);
  }

  // Renders this cell as image that fits in a square of the given side length
  // and as if it were the given distance away from some starting point with
  // respect to some maximum distance from the given point
  public WorldImage toImage(int sideLength, int distanceFromLoc, int maxDistanceFromLoc) {
    return this.toImage(sideLength,
        this.state.interpolatedColor(distanceFromLoc, maxDistanceFromLoc),
        Constants.INTERPOLATE_WALL_COLOR);
  }

  // returns a list of cells that are connected to this cell by open edges
  public ArrayList<ICell> adjacentOutCells() {

    ArrayList<CellEdge> allOutEdges = new ArrayList<CellEdge>(
        List.of(this.top, this.left, this.right, this.bottom));

    return new ArrayUtils().map(new ArrayUtils().filter(allOutEdges, new InTree()),
        new FirstConnectionThatIsnt(this));

  }

  // returns a list of cells that are connected to this cell
  public ArrayList<ICell> allAdjacentCells() {

    ArrayList<CellEdge> allOutEdges = new ArrayList<CellEdge>(
        List.of(this.top, this.left, this.right, this.bottom));

    return new ArrayUtils().map(allOutEdges, new FirstConnectionThatIsnt(this));
  }

  // EFFECT: changes the state of a cell, reassigning it to the given state
  public void enter(ICellState state) {
    this.state = state;
  }

  // returns the cell above this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeTop() {

    if (this.top.inTree()) {
      return this.top.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell below this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeBottom() {

    if (this.bottom.inTree()) {
      return this.bottom.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell left this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeLeft() {

    if (this.left.inTree()) {
      return this.left.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell right this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeRight() {

    if (this.right.inTree()) {
      return this.right.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopRight() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomRight() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomLeft() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopLeft() {
    return this;
  }
}

// Represents a hexagon cell in a maze. The hexagon itself
// is rotated 30 degrees with respect to the +X axis in the standard
// Cartesian coordinate system. That is, it is treated as a diamond
class HexCell implements ICell {

  private CellEdge topRight;
  private CellEdge topLeft;
  private CellEdge left;
  private CellEdge bottomLeft;
  private CellEdge bottomRight;
  private CellEdge right;
  private ICellState state;

  // Constructs a HexCell with connections to six DummyCells
  HexCell() {
    this.topRight = new CellEdge(this, new DummyCell());
    this.topLeft = new CellEdge(this, new DummyCell());
    this.left = new CellEdge(this, new DummyCell());
    this.bottomLeft = new CellEdge(this, new DummyCell());
    this.bottomRight = new CellEdge(this, new DummyCell());
    this.right = new CellEdge(this, new DummyCell());
    this.state = new UnseenState();
  }

  // Returns the new edge created connecting this HexCell to the given HexCell
  // EFFECT: Connects this cell to the given cell at its bottom right, and
  // connects
  // the given cell to this cell at its top left with the given weight
  CellEdge connectAtBottomRight(HexCell other, int weight) {
    CellEdge connection = new CellEdge(other, this, weight);
    this.bottomRight = connection;
    other.topLeft = connection;
    return connection;
  }

  // Returns the new edge created connecting this HexCell to the given HexCell
  // EFFECT: Connects this cell to the given cell at its bottom left, and connects
  // the given cell to this cell at its top right with the given weight
  CellEdge connectAtBottomLeft(HexCell other, int weight) {
    CellEdge connection = new CellEdge(other, this, weight);
    this.bottomLeft = connection;
    other.topRight = connection;
    return connection;
  }

  // Returns the new edge created connecting this HexCell to the given HexCell
  // EFFECT: Connects this cell to the given cell at its right, and connects
  // the given cell to this cell at its left with the given weight
  CellEdge connectAtRight(HexCell other, int weight) {
    CellEdge connection = new CellEdge(other, this, weight);
    this.right = connection;
    other.left = connection;
    return connection;
  }

  // returns a list of cells that are connected to this cell by open edges
  public ArrayList<ICell> adjacentOutCells() {

    ArrayList<CellEdge> allOutEdges = new ArrayList<CellEdge>(List.of(this.topLeft, this.topRight,
        this.left, this.right, this.bottomLeft, this.bottomRight));

    return new ArrayUtils().map(new ArrayUtils().filter(allOutEdges, new InTree()),
        new FirstConnectionThatIsnt(this));

  }

  // returns a list of cells that are connected to this cell
  public ArrayList<ICell> allAdjacentCells() {

    ArrayList<CellEdge> allOutEdges = new ArrayList<CellEdge>(List.of(this.topLeft, this.topRight,
        this.left, this.right, this.bottomLeft, this.bottomRight));

    return new ArrayUtils().map(allOutEdges, new FirstConnectionThatIsnt(this));
  }

  // EFFECT: changes the state of a cell, reassigning it to the given state
  public void enter(ICellState state) {
    this.state = state;
  }

  // Renders this Cell with the given wall color and background color
  private WorldImage toImageHelp(int sideLength, Color cellColor, Color wallColor) {
    WorldImage drawnCell = new HexagonImage(sideLength, OutlineMode.SOLID, cellColor);

    // The initial image is an unrotated hexagon. We rotate the image
    // counterclockwise from above
    // at the end after adding all of the walls, as the math is slightly
    // easier to work with. That being said, we must be careful which parts of the
    // unrotated
    // hexagon map to `this.topRight`, `this.topLeft`, etc.
    //
    // Currently, our `this.right` is at the top right of the unrotated hexagon
    // image and will become the right edge after rotation. The rest follows from
    // that

    if (!this.right.inTree()) {
      drawnCell = this.placeLine(0, sideLength, drawnCell);
    }

    if (!this.topRight.inTree()) {
      drawnCell = this.placeLine(Math.PI / 3, sideLength, drawnCell);
    }

    if (!this.topLeft.inTree()) {
      drawnCell = this.placeLine(2 * Math.PI / 3, sideLength, drawnCell);
    }

    if (!this.left.inTree()) {
      drawnCell = this.placeLine(Math.PI, sideLength, drawnCell);
    }

    if (!this.bottomLeft.inTree()) {
      drawnCell = this.placeLine(4 * Math.PI / 3, sideLength, drawnCell);
    }

    if (!this.bottomRight.inTree()) {
      drawnCell = this.placeLine(5 * Math.PI / 3, sideLength, drawnCell);
    }

    // RECALL: 30 degree towards the +Y axis from the +X axis. THIS IS IN GRAPHICS
    // SPACE (not absolute space), so it rotates CLOCKWISE with respect to absolute
    // space
    return new RotateImage(drawnCell, 30);
  }

  // Given an image to place a line into, produces a new image with
  // a new line segment perpendicular to the line described in polar coordinates
  // as `theta = baseTheta` of the given length at a distance `sqrt(3) / 2 *
  // sideLength`
  // away from the origin of the overlayed images (where their pinholes overlap).
  // When coupled with an image of a hexagon
  // of side length `sideLength`, the line segment produced will be snugly placed
  // at the border of the hexagon connecting the vertices of the hexagon
  // at `baseTheta` and `baseTheta + Math.PI / 3`
  private WorldImage placeLine(double baseTheta, int sideLength, WorldImage hexImage) {

    // The way this works is as follows:
    // There are two coordinate systems at play here: the standard
    // Cartesian coordinate system (SCS), with +X pointing right and +Y pointing up,
    // and
    // the standard coordinate system for 2D computer graphics (+X pointing right
    // and +Y pointing down) (SGS). The latter is used because it is easier to
    // conceptualize.
    //
    // To draw a line in the correct spot in a hexagon, first note
    // that we start by imagining the hexagon in the SCS. `baseTheta` is then
    // the angle with respect to the positive +X axis that one end of the line
    // will lie at. The other end will be at angle `baseTheta + Math.PI / 3`, or
    // 60 degrees ccw further than `baseTheta`.
    //
    // To get the line to point at `baseTheta + Math.PI / 3`, we note
    // first that the line constructor is made to point in its own SGS.
    // Hence, we must convert to that system from the absolute SCS.
    //
    // Notice that the side of an unrotated hexagon is parallel to the side
    // 120 degrees ccw to it (and also 240, but 120 suffices here); hence, if
    // we point the line at `baseTheta + 2 * Math.PI / 3` FROM THE ORIGIN (since
    // that
    // is how the line is constructed according to the docs, converting from SCS
    // to its SGS), the line will now point in the same direction `baseTheta +
    // Math.PI / 6`
    // as if it were pointing at that point from `baseTheta`.
    //
    // Finally, because the Image library uses Posns, a conversion is made
    // with as little error as possible (we don't use `Posn`s throughout because
    // of floating-point inaccuracy accrual).

    I2DCoordinateSystem whereToPlace = new RotatedCoordinateSystem(baseTheta + Math.PI / 6,
        new StandardCoordinateSystem());
    I2DCoordinateSystem whereToDirect = new RotatedCoordinateSystem(baseTheta + 2 * Math.PI / 3,
        new StandardCoordinateSystem());
    I2DCoordinateSystem graphics = new GraphicsCoordinateSystem();

    Point2D locToDirectAbsolute = whereToDirect.convertToWorld(new Point2D(sideLength, 0));
    Point2D lineImageToPlace = graphics.convertToThis(locToDirectAbsolute);
    Posn approxLocationToPlace = new Posn((int) (Math.round(lineImageToPlace.x)),
        (int) (Math.round(lineImageToPlace.y)));

    WorldImage line = new LineImage(approxLocationToPlace, Color.BLACK);

    // -loc.x because we need to move it back, loc.y to move it UP (in graphics
    // system)

    Point2D locToPlace = whereToPlace.convertToWorld(new Point2D(sideLength * Math.sqrt(3) / 2, 0));

    return new OverlayOffsetImage(line, -locToPlace.x, locToPlace.y, hexImage);
  }

  // renders this cell as an image
  public WorldImage toImage(int sideLength) {
    return this.toImageHelp(sideLength, this.state.standardColor(), Constants.WALL_COLOR);
  }

  // Renders this cell as image that fits in a square of the given side length
  // and as if it were the given distance away from some starting point with
  // respect to some maximum distance from the given point
  public WorldImage toImage(int sideLength, int distanceFromLoc, int maxDistanceFromLoc) {
    return this.toImageHelp(sideLength,
        this.state.interpolatedColor(distanceFromLoc, maxDistanceFromLoc),
        Constants.INTERPOLATE_WALL_COLOR);
  }

  // returns the cell above this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeTop() {
    return this;
  }

  // returns the cell below this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeBottom() {
    return this;
  }

  // returns the cell left this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeLeft() {
    if (this.left.inTree()) {
      return this.left.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell right this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeRight() {
    if (this.right.inTree()) {
      return this.right.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopRight() {
    if (this.topRight.inTree()) {
      return this.topRight.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomRight() {
    if (this.topRight.inTree()) {
      return this.topRight.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomLeft() {
    if (this.topRight.inTree()) {
      return this.topRight.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopLeft() {
    if (this.topRight.inTree()) {
      return this.topRight.firstConnectionThatIsnt(this);
    } else {
      return this;
    }
  }
}

//represents a placeholder cell that can either hold the place of a future
//connection or form the boundary of a maze
class DummyCell implements ICell {

  // renders the cell as an image
  public WorldImage toImage(int sideLength) {
    return new EmptyImage();
  }

  // Renders this cell as image that fits in a square of the given side length
  // and as if it were the given distance away from some starting point with
  // respect to some maximum distance from the given point
  public WorldImage toImage(int sideLength, int distanceFromLoc, int maxDistanceFromLoc) {
    return new EmptyImage();
  }

  // returns a list of cells that are connected to this cell by open edges
  public ArrayList<ICell> adjacentOutCells() {
    return new ArrayList<ICell>();
  }

  // returns a list of cells that are connected to this cell
  public ArrayList<ICell> allAdjacentCells() {
    return new ArrayList<ICell>();
  }

  // EFFECT: changes the state of a cell, reassigning it to the given state
  // Has no effect on DummyCells, as they lack states, so an exception is thrown
  public void enter(ICellState state) {
    throw new UnsupportedOperationException(
        "DummyCells lack states, so their" + " states cannot be modified.");
  }

  // returns the cell above this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeTop() {
    return this;
  }

  // returns the cell below this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeBottom() {
    return this;
  }

  // returns the cell left this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeLeft() {
    return this;
  }

  // returns the cell right this cell, if there is an open connection to it,
  // returning this cell otherwise
  public ICell maybeRight() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopRight() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomRight() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeBottomLeft() {
    return this;
  }

  // returns the cell to the top right of this cell, if there is an open
  // connection to it,
  // returning this cell otherwise
  public ICell maybeTopLeft() {
    return this;
  }

}

// represents an undirected connection between two ICells in a maze
final class CellEdge implements Comparable<CellEdge> {
  public final ICell first;
  public final ICell second;
  public final int weight;
  private IEdgeState state;

  // Construct a new CellEdge representing a connection
  // between two ICells of the given weight
  CellEdge(ICell first, ICell second, int weight) {
    if (weight < 0) {
      throw new IllegalArgumentException(
          "Cannot represent a connection between two ICells with a negative weight");
    }
    this.first = first;
    this.second = second;
    this.weight = weight;
    this.state = new ClosedState();
  }

  // constructs a CellEdge with default weight 0
  CellEdge(ICell first, ICell second) {
    this(first, second, 0);
  }

  // is this CellEdge a member of the minimum spanning tree of all open paths?
  boolean inTree() {
    return this.state.isOpenPath();
  }

  // EFFECT: Erases knowledge from this edge that it is
  // part of the minimum spanning tree through the cells
  // of the Maze of which our two cells are apart of
  void removeFromTree() {
    this.state = new ClosedState();
  }

  // EFFECT: Provides knowledge to this edge that it is
  // now part of the minimum spanning tree through the cells
  // of the Maze of which our two cells are apart of
  void becomeEdgeInTree() {
    this.state = new OpenState();
  }

  // Compares this edge with the other edge by
  // comparing the weight of this edge with the given one
  public int compareTo(CellEdge other) {
    return this.weight - other.weight;
  }

  // Returns the first connection to this edge that isn't the given cell
  // using intentional equality. If the edge is such that it points to the
  // (intentionally) same cell, an exception is thrown
  public ICell firstConnectionThatIsnt(ICell thatCell) {
    if (this.first == thatCell) {
      if (this.second == thatCell) {
        throw new IllegalArgumentException(
            "This CellEdge points to the given " + "cell in both directions.");
      } else {
        return this.second;
      }
    } else {
      return this.first;
    }
  }
}

// Represents a state that a cell can have
interface ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  Color standardColor();

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  Color interpolatedColor(int lower, int upper);
}

//represents the state of a cell that has not yet been seen in the search
final class UnseenState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_UNSEEN_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    double ratio = (double) lower / (double) upper;
    return new ColorUtils().interpolate(ratio, Constants.CELL_CLOSE_COLOR,
        Constants.CELL_FAR_COLOR);
  }
}

//represents the state of a cell that has been seen in the search but is
//not a part of the completed solution
final class SeenState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_SEEN_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    return this.standardColor();
  }
}

//represents the state of a cell that is in the completed solution
final class PathState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_PATH_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    return this.standardColor();
  }
}

//represents the state of a cell that is the starting point of the maze
final class StartState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_START_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    double ratio = (double) lower / (double) upper;
    return new ColorUtils().interpolate(ratio, Constants.CELL_CLOSE_COLOR,
        Constants.CELL_FAR_COLOR);
  }
}

//represents the state of a cell that is the ending point of the maze
final class FinishState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_FINISH_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    double ratio = (double) lower / (double) upper;
    return new ColorUtils().interpolate(ratio, Constants.CELL_CLOSE_COLOR,
        Constants.CELL_FAR_COLOR);
  }
}

//represents the state of a cell that is the our current "head", meaning
//that it is the cell most recently seen
final class HeadState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_HEAD_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    return this.standardColor();
  }
}

//represents the state of a cell that is not our current head, but is instead
//in the collection of cells to be searched (only used in breadth first search
//animation)
final class FutureHeadState implements ICellState {
  // produces the color associated with this state when no color-affecting
  // alterations are applied
  public Color standardColor() {
    return Constants.CELL_FUTUREHEAD_COLOR;
  }

  // Produces the color associated with this state when we are coloring
  // based on interpolating between colors based on distance
  public Color interpolatedColor(int lower, int upper) {
    return this.standardColor();
  }
}

//represents a state that a CellEdge can have
interface IEdgeState {
  // is this a state of a CellEdge that describes that Edge as a part of
  // the minimum spanning tree?
  boolean isOpenPath();
}

//represents the state of a CellEdge that represents a wall in the maze
final class ClosedState implements IEdgeState {
  // is this a state of a CellEdge that describes that Edge as a part of
  // the minimum spanning tree?
  public boolean isOpenPath() {
    return false;
  }
}

//represents the state of a CellEdge that represents an open path in the maze
final class OpenState implements IEdgeState {
  // is this a state of a CellEdge that describes that Edge as a part of
  // the minimum spanning tree?
  public boolean isOpenPath() {
    return true;
  }
}

//represents an object that can create horizontal and vertical edge weights
// according to two different distributions
final class EdgeWeightGen {

  private final Random hRand;
  private final Random vRand;
  private final int hMax;
  private final int vMax;

  // constructs an Edge Weight Generator with maximum weights for vertical
  // and horizontal weights as well as random generator seeds for each
  // weight distribution
  EdgeWeightGen(int hMax, int vMax, int hSeed, int vSeed) {
    if (hMax < 0 || vMax < 0) {
      throw new IllegalArgumentException(
          "Maximum weights must be greater than" + " or equal to zero.");
    }
    this.hRand = new Random(hSeed);
    this.vRand = new Random(vSeed);
    this.hMax = hMax;
    this.vMax = vMax;
  }

  // Constructs a EdgeWeightGen with specific seeds for non-variable
  // random seeds
  EdgeWeightGen(int hMax, int vMax) {
    this(hMax, vMax, 20, 20);
  }

  // produces the next horizontal weight
  int nextHorizontalWeight() {
    return this.hRand.nextInt(hMax + 1);
  }

  // produces the next vertical weight
  int nextVerticalWeight() {
    return this.vRand.nextInt(vMax + 1);
  }
}

//represents utility functions on ArrayLists
class ArrayUtils {

  // Computes a single value from the given list by combining each item in the
  // list with the accumulated result so far, moving from left to right, starting
  // with the base value. If the list is empty, the base value is returned.
  <T, U> U foldl(ArrayList<T> list, U base, BiFunction<T, U, U> func) {
    U currResult = base;

    for (T element : list) {
      currResult = func.apply(element, currResult);
    }

    return currResult;
  }

  // returns the first item in the list that satisfies the predicate,
  // throwing an exception if none satisfy
  <T> T firstToSatisfy(ArrayList<T> list, Predicate<T> predicate) {
    if (list.isEmpty()) {
      throw new IllegalArgumentException("Cannot find an item in an empty list");
    }

    for (T element : list) {
      if (predicate.test(element)) {
        return element;
      }
    }

    throw new IllegalArgumentException(
        "No argument in the list satisfied the given" + " function.");
  }

  // returns an ArrayList of all items in the list that satisfy the given
  // predicate
  <T> ArrayList<T> filter(ArrayList<T> list, Predicate<T> predicate) {
    ArrayList<T> soFar = new ArrayList<T>();
    for (T t : list) {
      if (predicate.test(t)) {
        soFar.add(t);
      }
    }
    return soFar;
  }

  // applies a given function to every item in an ArrayList, returning an
  // ArrayList of each result
  <T, U> ArrayList<U> map(ArrayList<T> list, Function<T, U> func) {
    ArrayList<U> soFar = new ArrayList<U>();
    for (T t : list) {
      soFar.add(func.apply(t));
    }
    return soFar;
  }
}

// A utilities class for doubles
class DoubleUtils {

  // Produces an integer between the two given values
  // or the value itself if it is between the two thresholds
  double clamp(double value, double lower, double higher) {
    if (lower > higher) {
      throw new IllegalArgumentException(
          "Cannot clamp a value in a range where the lower value is greater than the higher one");
    }
    if (value < lower) {
      return lower;
    } else if (value > higher) {
      return higher;
    } else {
      return value;
    }
  }
}

// A utilities class for colors
class ColorUtils {

  // Produces a color that is the linear interpolation of the two given colors
  // given the relative weight between 0 and 1
  Color interpolate(double weight, Color col1, Color col2) {
    double col1Red = col1.getRed();
    double col1Green = col1.getGreen();
    double col1Blue = col1.getBlue();

    double col2Red = col2.getRed();
    double col2Green = col2.getGreen();
    double col2Blue = col2.getBlue();

    double clampedWeight = new DoubleUtils().clamp(weight, 0.0, 1.0);

    double newRed = new DoubleUtils().clamp(
        (1.0 - clampedWeight) * col1Red + clampedWeight * col2Red, Math.min(col1Red, col2Red),
        Math.max(col1Red, col2Red));
    double newGreen = new DoubleUtils().clamp(
        (1.0 - clampedWeight) * col1Green + clampedWeight * col2Green,
        Math.min(col1Green, col2Green), Math.max(col1Green, col2Green));
    double newBlue = new DoubleUtils().clamp(
        (1.0 - clampedWeight) * col1Blue + clampedWeight * col2Blue, Math.min(col1Blue, col2Blue),
        Math.max(col1Blue, col2Blue));

    return new Color((int) newRed, (int) newGreen, (int) newBlue);
  }
}

//a predicate to determine if an item is in an addFind
class AddFindContains<T> implements Predicate<T> {

  private final AddFind<T> addFind;

  AddFindContains(AddFind<T> addFind) {
    this.addFind = addFind;
  }

  // is the given element in our addFind?
  public boolean test(T t) {
    return this.addFind.contains(t);
  }

}

//a predicate to negate another predicate
class NotPred<T> implements Predicate<T> {

  private final Predicate<T> pred;

  NotPred(Predicate<T> pred) {
    this.pred = pred;
  }

  public boolean test(T t) {
    return !this.pred.test(t);
  }

}

//is the given cellEdge open and in the minimum spanning tree?
class InTree implements Predicate<CellEdge> {

  public boolean test(CellEdge edge) {
    return edge.inTree();
  }

}

//returns the first connection to a CellEdge that isn't the given cell
class FirstConnectionThatIsnt implements Function<CellEdge, ICell> {
  ICell thisEdge;

  FirstConnectionThatIsnt(ICell thisEdge) {
    this.thisEdge = thisEdge;
  }

  public ICell apply(CellEdge edge) {
    return edge.firstConnectionThatIsnt(thisEdge);
  }
}

// Describes the custom maze world the user interacts with and has fun with
final class MazeWorld extends World {
  private Maze worldMaze;
  private ASearchSnapshot searchSnapshot;
  private boolean isPaused;
  // are we running DFS? (or are we running BFS)
  // ignored if isManual is true
  private boolean isDepthSearch;
  // is the user manually solving the maze?
  private boolean isManual;
  // is the maze hexagonal?
  private boolean isHex;
  // have we found our intended finishing cell?
  private boolean foundTheEnd;
  // a random number generator for building mazes
  // by seed value
  private Random seedGen;
  // the number of rows in the maze
  private final int cellsPerRow;
  // the number of columns in the maze
  private final int cellsPerColumn;

  // Are we rendering the scene with distance from the start?
  private boolean renderAsDistanceFromStart;

  // Are we rendering the scene with distance from the end?
  private boolean renderAsDistanceFromFinish;

  // Are we showing the seen paths
  private boolean showSeen;

  // The width of the window into which
  // the world is drawn
  private final int windowWidth;

  // The width of the window into which
  // the world is drawn
  private final int windowHeight;

  // construct a new maze with `mazeWidth` number of cells per
  // row and `mazeHeight` number of cells per column and
  // a special seed that determines the order in which the user
  // encounters new mazes
  MazeWorld(int cellsPerRow, int cellsPerColumn, int windowWidth, int windowHeight,
      int seedGenerationSeed, boolean isHex) {

    // Ensure that there is at least 1 cell in each row and column
    if (cellsPerRow <= 0) {
      throw new IllegalArgumentException(
          "A maze must have at least one cell per row. Given " + Integer.toString(cellsPerRow));
    }

    if (cellsPerColumn <= 0) {
      throw new IllegalArgumentException("A maze must have at least one cell per column. Given "
          + Integer.toString(cellsPerColumn));
    }

    // Ensure that the window has a nonnegaitve width and height
    if (windowWidth <= 0) {
      throw new IllegalArgumentException("A maze must have a positive window width");
    }

    if (windowHeight <= 0) {
      throw new IllegalArgumentException("A maze must have a positive window height");
    }

    this.seedGen = new Random(seedGenerationSeed);
    this.isPaused = false;
    this.isDepthSearch = true;
    this.isManual = false;
    this.isHex = isHex;
    this.foundTheEnd = false;
    this.cellsPerColumn = cellsPerColumn;
    this.cellsPerRow = cellsPerRow;
    this.renderAsDistanceFromStart = false;
    this.showSeen = true;
    this.windowWidth = windowWidth;
    this.windowHeight = windowHeight;

    this.worldMaze = new Maze(cellsPerRow, cellsPerColumn, this.seedGen.nextInt(), this.isHex);
    this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
  }

  // Produces the scene that is displayed by `bigBang` to the user
  public WorldScene makeScene() {
    WorldScene scene = new WorldScene(this.windowWidth, this.windowHeight);

    if (this.renderAsDistanceFromStart && !this.renderAsDistanceFromFinish) {
      this.worldMaze.drawInSceneInterpolated(scene, this.windowWidth, this.windowHeight, true);
    } else if (!this.renderAsDistanceFromStart && this.renderAsDistanceFromFinish) {
      this.worldMaze.drawInSceneInterpolated(scene, this.windowWidth, this.windowHeight, false);
    } else if (this.renderAsDistanceFromStart && this.renderAsDistanceFromFinish) {
      throw new IllegalStateException("A maze cannot show it's distances from the start and "
          + "distances from the end at the same time.");
    } else {
      this.worldMaze.drawInScene(scene, this.windowWidth, this.windowHeight);
    }

    return scene;
  }

  // EFFECT: modifies the MazeWorld based on the state of the Maze and
  // SearchSnapshot,
  // either breaking down one wall to create the maze or taking one step into a
  // new cell
  // to solve the maze
  public void onTick() {

    if (this.foundTheEnd) {
      return;
    }

    if (this.isPaused) {
      return;
    }

    if (!this.worldMaze.doneBreakingWalls()) {
      this.worldMaze.breakOneWall();
      return;
    }

    if (this.searchSnapshot.hasMoreToSee()) {

      this.foundTheEnd = this.searchSnapshot.searchOneStep();
    }

    if (this.showSeen) {
      this.searchSnapshot.makeAllSeen();
    } else {
      this.searchSnapshot.makeOnlyHead();
    }

    // if we've found the end, we'll show the final path
    if (this.foundTheEnd) {
      this.searchSnapshot.makeAllInPath();
    }

  }

  // EFFECT: Performs one of numerous actions depending
  // on what key the user presses (see the gameplay manual
  // on what each key does), from making a new maze
  // to starting a new search through the maze
  public void onKeyEvent(String key) {

    if (key.equals(" ")) {
      this.togglePause();
    } else if (key.equals("b")) {
      this.setBreadthFirst();
    } else if (key.equals("p")) {
      this.setDepthFirst();
    } else if (key.equals("n")) {
      this.newMaze();
    } else if (key.equals("r")) {
      this.newRowMaze();
    } else if (key.equals("c")) {
      this.newColumnMaze();
    } else if (key.equals("!")) {
      // Skip the maze-wall knocking down
      this.worldMaze.breakAllWalls();
    } else if (key.equals("k") && this.worldMaze.doneBreakingWalls()) {
      this.toggleDistanceFromStart();
    } else if (key.equals("l") && this.worldMaze.doneBreakingWalls()) {
      this.toggleDistanceFromFinish();
    } else if (key.equals("m")) {
      this.setManual();
    } else if (key.equals("up") && this.isManual && !isHex) {
      this.searchSnapshot.maybeMoveUp();
    } else if (key.equals("down") && this.isManual && !isHex) {
      this.searchSnapshot.maybeMoveDown();
    } else if (key.equals("left") && this.isManual && !isHex) {
      this.searchSnapshot.maybeMoveLeft();
    } else if (key.equals("right") && this.isManual && !isHex) {
      this.searchSnapshot.maybeMoveRight();
    } else if (key.equals("-")) {
      this.showSeen = !this.showSeen;
    } else if (key.equals("w") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveTopLeft();
    } else if (key.equals("e") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveTopRight();
    } else if (key.equals("a") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveLeft();
    } else if (key.equals("d") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveRight();
    } else if (key.equals("z") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveBottomLeft();
    } else if (key.equals("x") && this.isHex && this.isManual) {
      this.searchSnapshot.maybeMoveBottomRight();
    } else if (key.equals("q")) {
      this.newHexMaze();
    }
  }

  // pauses or unpauses the mazeWorld
  private void togglePause() {

    this.isPaused = !this.isPaused;
  }

  // begins searching the maze with a new Breadth First Search
  private void setBreadthFirst() {
    this.isDepthSearch = false;
    this.isManual = false;
    this.searchSnapshot.makeAllUnseen();
    this.searchSnapshot = this.worldMaze.newBreadthFirstSnapshot();
    this.foundTheEnd = false;
  }

  // begins searching the maze with a new Depth First Search
  private void setDepthFirst() {
    this.isDepthSearch = true;
    this.isManual = false;
    this.searchSnapshot.makeAllUnseen();
    this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
    this.foundTheEnd = false;
  }

  // creates a new maze without bias
  private void newMaze() {

    this.renderAsDistanceFromStart = false;
    this.renderAsDistanceFromFinish = false;

    this.worldMaze = new Maze(this.cellsPerRow, this.cellsPerColumn, this.seedGen.nextInt(),
        this.isHex);

    if (this.isManual) {
      this.searchSnapshot = this.worldMaze.newManualSnapshot();
    } else {
      if (this.isDepthSearch) {
        this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
      } else {
        this.searchSnapshot = this.worldMaze.newBreadthFirstSnapshot();
      }
    }

    this.foundTheEnd = false;
  }

  // creates a new maze, toggling hexagonal mazes
  private void newHexMaze() {

    this.isHex = !this.isHex;

    this.renderAsDistanceFromStart = false;
    this.renderAsDistanceFromFinish = false;

    this.worldMaze = new Maze(this.cellsPerRow, this.cellsPerColumn, this.seedGen.nextInt(),
        this.isHex);

    if (this.isManual) {
      this.searchSnapshot = this.worldMaze.newManualSnapshot();
    } else {
      if (this.isDepthSearch) {
        this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
      } else {
        this.searchSnapshot = this.worldMaze.newBreadthFirstSnapshot();
      }
    }

    this.foundTheEnd = false;

  }

  // creates a new maze with horizontal bias
  private void newRowMaze() {

    this.renderAsDistanceFromStart = false;
    this.renderAsDistanceFromFinish = false;

    this.worldMaze = new Maze(this.cellsPerRow, this.cellsPerColumn, this.seedGen.nextInt(), true,
        false, this.isHex);

    if (this.isManual) {
      this.searchSnapshot = this.worldMaze.newManualSnapshot();
    } else {
      if (this.isDepthSearch) {
        this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
      } else {
        this.searchSnapshot = this.worldMaze.newBreadthFirstSnapshot();
      }
    }

    this.foundTheEnd = false;

  }

  // creates a new maze with vertical bias
  private void newColumnMaze() {

    this.renderAsDistanceFromStart = false;
    this.renderAsDistanceFromFinish = false;

    this.worldMaze = new Maze(this.cellsPerRow, this.cellsPerColumn, this.seedGen.nextInt(), false,
        true, this.isHex);

    if (this.isManual) {
      this.searchSnapshot = this.worldMaze.newManualSnapshot();
    } else {
      if (this.isDepthSearch) {
        this.searchSnapshot = this.worldMaze.newDepthFirstSnapshot();
      } else {
        this.searchSnapshot = this.worldMaze.newBreadthFirstSnapshot();
      }
    }

    this.foundTheEnd = false;
  }

  // shows or hides the distances from the start of the maze
  private void toggleDistanceFromStart() {
    this.renderAsDistanceFromStart = !this.renderAsDistanceFromStart;
    this.renderAsDistanceFromFinish = false;
  }

  // shows or hides the distances from the start of the maze
  private void toggleDistanceFromFinish() {
    this.renderAsDistanceFromFinish = !this.renderAsDistanceFromFinish;
    this.renderAsDistanceFromStart = false;
  }

  // begins searching the maze with a manual search
  private void setManual() {
    this.isManual = true;
    this.searchSnapshot.makeAllUnseen();
    this.searchSnapshot = this.worldMaze.newManualSnapshot();
    this.foundTheEnd = false;
  }
}

// Represents a double-precision point in space.
// This is effectively a value type, and as such openly has its
// fields available to all users of the class
final class Point2D {
  public final double x;
  public final double y;

  // Construct a new point at the given X and Y values
  Point2D(double x, double y) {
    this.x = x;
    this.y = y;
  }

  // Construct a new Point2D from a Posn
  Point2D(Posn pos) {
    this(pos.x, pos.y);
  }

  // Is the given point within a box of the
  // given side length if it were centered at this point?
  public boolean wouldContainInBox(double sideLength, Point2D other) {
    return Math.abs(other.x - this.x) < sideLength && Math.abs(other.y - this.y) < sideLength;
  }
}

// Represents a 2D coordinate system in space fixed at the origin
// that is sufficient for converting `Posn`s in space
interface I2DCoordinateSystem {

  // Given a `Point2D` that describes a point in the absolute
  // coordinate system (i.e. in the standard +X pointing right and the
  // +Y pointing up each of unit scale), computes the coordinates of that
  // point in this coordinate system
  Point2D convertToThis(Point2D worldPos);

  // Converts the given Point2D, specified as coordinates in this
  // coordinate system, to coordinates in the absolute frame
  // (i.e. in the standard +X pointing right and the +Y pointing up each of unit
  // scale)
  Point2D convertToWorld(Point2D pointInSystem);
}

// Represents the standard 2D coordinate system seen typically
// in math classes (+X pointing right, and +Y pointing up)
class StandardCoordinateSystem implements I2DCoordinateSystem {
  // Given a `Point2D` that describes a point in the absolute
  // coordinate system (i.e. in the standard +X pointing right and the
  // +Y pointing up each of unit scale), computes the coordinates of that
  // point in this coordinate system
  public Point2D convertToThis(Point2D worldPos) {
    return worldPos;
  }

  // Converts the given Point2D, specified as coordinates in this
  // coordinate system, to coordinates in the absolute frame
  // (i.e. in the standard +X pointing right and the +Y pointing up each of unit
  // scale)
  public Point2D convertToWorld(Point2D pointInSystem) {
    return pointInSystem;
  }
}

// Represents the standard 2D coordinate system seen typically
// in graphics rendering (+X pointing right, and +Y pointing down)
class GraphicsCoordinateSystem implements I2DCoordinateSystem {

  // Given a `Point2D` that describes a point in the absolute
  // coordinate system (i.e. in the standard +X pointing right and the
  // +Y pointing up each of unit scale), computes the coordinates of that
  // point in this coordinate system
  public Point2D convertToThis(Point2D worldPos) {
    return new Point2D(worldPos.x, -worldPos.y);
  }

  // Converts the given Point2D, specified as coordinates in this
  // coordinate system, to coordinates in the absolute frame
  // (i.e. in the standard +X pointing right and the +Y pointing up each of unit
  // scale)
  public Point2D convertToWorld(Point2D pointInSystem) {
    return new Point2D(pointInSystem.x, -pointInSystem.y);
  }
}

// Represents a fixed, rotated version of a 2D coordinate system
class RotatedCoordinateSystem implements I2DCoordinateSystem {

  // The counterclockwise angle, when viewed from above measured in radians, that
  // the coordinate
  // system wrapped by this system is rotated at
  private final double theta;

  // The coordinate system that undergoes this rotation
  private final I2DCoordinateSystem system;

  // Create a RotatedCoordinateSystem by rotating the given
  // coordinate system by the given value
  RotatedCoordinateSystem(double theta, I2DCoordinateSystem system) {
    this.theta = theta;
    this.system = system;
  }

  // Given a `Point2D` that describes a point in the absolute
  // coordinate system (i.e. in the standard +X pointing right and the
  // +Y pointing up each of unit scale), computes the coordinates of that
  // point in this coordinate system
  public Point2D convertToThis(Point2D worldPos) {

    // First, convert the point into the system
    // that is rotated
    Point2D inSystem = this.system.convertToThis(worldPos);

    // Now, `inSystem` is effectively a point in the
    // standard coordinate system from the perspective of
    // `this.system`. All that remains it to convert
    // that point into this system. Since we are a coordinate
    // system rotated with respect to `this.system`, we must
    // rotate the point by `-this.theta` in order to convert it to
    // this coordinate system because a change of basis by a matrix
    // `A` whose column vectors are basis vectors has an effect of
    // `inverse A` to get the coordinates
    return this.rotated(inSystem, -this.theta);
  }

  // Produces a new Point2D by rotating through by the given
  // angle with respect to the +X axis, assuming that the point
  // is viewed from a standard frame of reference
  private Point2D rotated(Point2D pos, double theta) {
    double xComp = pos.x * Math.cos(theta) - pos.y * Math.sin(theta);
    double yComp = pos.x * Math.sin(theta) + pos.y * Math.cos(theta);
    return new Point2D(xComp, yComp);
  }

  // Converts the given Point2D, specified as coordinates in this
  // coordinate system, to coordinates in the absolute frame
  // (i.e. in the standard +X pointing right and the +Y pointing up each of unit
  // scale)
  public Point2D convertToWorld(Point2D pointInSystem) {
    // First, we convert to the frame of the
    // system that this system represents are rotated
    Point2D posToUnderlyingSystem = this.rotated(pointInSystem, this.theta);

    // And simply convert that point to world space
    return this.system.convertToWorld(posToUnderlyingSystem);
  }
}

//represents a function object appending lists together
final class AppendList<T> implements BiFunction<ArrayList<T>, ArrayList<T>, ArrayList<T>> {
  public ArrayList<T> apply(ArrayList<T> next, ArrayList<T> prev) {
    // Create a copy of the contents so as to not modify the lists
    ArrayList<T> prevCopy = new ArrayList<T>(prev);
    prevCopy.addAll(next);
    return prevCopy;
  }
}

// Sorts CellEdges by their weights
final class ByEdgeWeight implements Comparator<CellEdge> {
  public int compare(CellEdge edge1, CellEdge edge2) {
    return edge1.compareTo(edge2);
  }
}

// FOR TESTING:

// Maps strings to their lengths
class ToLength implements Function<String, Integer> {
  public Integer apply(String s) {
    return s.length();
  }
}

// Maps inters to their additive inverses
class IntegerAdditiveInverse implements Function<Integer, Integer> {
  public Integer apply(Integer value) {
    return -value;
  }
}

// Maps a two integers to their combined value by adding them together
class TotalValue implements BiFunction<Integer, Integer, Integer> {
  public Integer apply(Integer value1, Integer value2) {
    return value1 + value1;
  }
}

// Combines two strings together
class CombineStrings implements BiFunction<String, String, String> {
  public String apply(String s1, String s2) {
    return s2 + s1;
  }
}

// Determines if an integer is positive
class IsPositive implements Predicate<Integer> {
  public boolean test(Integer value) {
    return value > 0;
  }
}

// Determines if an integer is less than or equal to 0
class IsNegativeOrZero implements Predicate<Integer> {
  public boolean test(Integer value) {
    return value <= 0;
  }
}

// Is the Posn in the first quadrant
class InFirstQuadrant implements Predicate<Posn> {
  public boolean test(Posn posn) {
    return posn.x > 0 && posn.y > 0;
  }
}

// Examples and tests for the maze game
class ExamplesMazes {

  // Images of cells of size 100x100 with all of its walls
  // in the various states
  WorldImage boxCellImageUnseen100x100;
  WorldImage boxCellImageSeen100x100;
  WorldImage boxCellImageStart100x100;
  WorldImage boxCellImageFinish100x100;
  WorldImage boxCellImagePath100x100;

  // Images of cells of size 100x100 with all of its walls
  // in the various states with interpolated color values
  WorldImage boxCellImageUnseenInterpolated100x100;
  WorldImage boxCellImageSeenInterpolated100x100;
  WorldImage boxCellImageStartInterpolated100x100;
  WorldImage boxCellImageFinishInterpolated100x100;
  WorldImage boxCellImagePathInterpolated100x100;

  // Images of cells of size 150x150 with all of its walls
  // in the various states
  WorldImage boxCellImageUnseen150x150;
  WorldImage boxCellImageSeen150x150;
  WorldImage boxCellImageStart150x150;
  WorldImage boxCellImageFinish150x150;
  WorldImage boxCellImagePath150x150;

  // Images of cells of size 100x100 with all of its walls
  // in the various states with interpolated color values
  WorldImage boxCellImageUnseenInterpolated150x150;
  WorldImage boxCellImageSeenInterpolated150x150;
  WorldImage boxCellImageStartInterpolated150x150;
  WorldImage boxCellImageFinishInterpolated150x150;
  WorldImage boxCellImagePathInterpolated150x150;

  // Images of hexagonal cells of size 100x100 with all of its walls
  // in the various states
  WorldImage hexCellImageUnseen100x100;
  WorldImage hexCellImageSeen100x100;
  WorldImage hexCellImageStart100x100;
  WorldImage hexCellImageFinish100x100;
  WorldImage hexCellImagePath100x100;

  // Images of hexagonal cells of size 100x100 with all of its walls
  // in the various states with interpolated color values
  WorldImage hexCellImageUnseenInterpolated100x100;
  WorldImage hexCellImageSeenInterpolated100x100;
  WorldImage hexCellImageStartInterpolated100x100;
  WorldImage hexCellImageFinishInterpolated100x100;
  WorldImage hexCellImagePathInterpolated100x100;

  // The maze drawn out looks like
  // __
  // |__|
  //
  Maze oneByOneMaze;

  // The maze drawn out looks like
  // __ __
  // | __|
  // |__ __|
  //
  Maze twoByTwoMaze;

  // The maze drawn out looks like
  // __ __ __
  // | |__ |
  // | _____|
  // |__ __ __|
  //
  Maze threeByThreeMaze;

  // The maze drawn out looks like
  //
  // HEX -- HEX HEX
  // / \
  // / \
  // HEX HEX -- HEX
  // \ / /
  // \ / /
  // HEX -- HEX HEX
  //
  Maze threeByThreeHexMaze;

  // Given an image to place a line into, produces a new image with
  // a new line segment perpendicular to the line described in polar coordinates
  // as `theta = baseTheta` of the given length at a distance `sqrt(3) / 2 *
  // sideLength`
  // away from the origin of the overlayed images. See `HexCell.toImage` for an
  // extended discussion
  private WorldImage placeLine(double baseTheta, int sideLength, WorldImage hexImage) {
    I2DCoordinateSystem whereToPlace = new RotatedCoordinateSystem(baseTheta + Math.PI / 6,
        new StandardCoordinateSystem());
    I2DCoordinateSystem whereToDirect = new RotatedCoordinateSystem(baseTheta + 2 * Math.PI / 3,
        new StandardCoordinateSystem());
    I2DCoordinateSystem graphics = new GraphicsCoordinateSystem();

    Point2D locToDirectAbsolute = whereToDirect.convertToWorld(new Point2D(sideLength, 0));
    Point2D lineImageToPlace = graphics.convertToThis(locToDirectAbsolute);
    Posn approxLocationToPlace = new Posn((int) (Math.round(lineImageToPlace.x)),
        (int) (Math.round(lineImageToPlace.y)));

    WorldImage line = new LineImage(approxLocationToPlace, Color.BLACK);

    // -loc.x because we need to move it back, loc.y to move it UP (in graphics
    // system)

    Point2D locToPlace = whereToPlace.convertToWorld(new Point2D(sideLength * Math.sqrt(3) / 2, 0));

    return new OverlayOffsetImage(line, -locToPlace.x, locToPlace.y, hexImage);
  }

  // Constructs new hexagon image of the given side length,
  // background color, wall color, and whether or not to
  // draw certain walls
  private WorldImage boxedHexagonImage(int sideLength, Color backgroundColor, Color wallColor,
      boolean topRight, boolean topLeft, boolean left, boolean right, boolean bottomLeft,
      boolean bottomRight) {

    WorldImage drawnCell = new HexagonImage(sideLength, OutlineMode.SOLID, backgroundColor);
    if (right) {
      drawnCell = this.placeLine(0, sideLength, drawnCell);
    }

    if (topRight) {
      drawnCell = this.placeLine(Math.PI / 3, sideLength, drawnCell);
    }

    if (topLeft) {
      drawnCell = this.placeLine(2 * Math.PI / 3, sideLength, drawnCell);
    }

    if (left) {
      drawnCell = this.placeLine(Math.PI, sideLength, drawnCell);
    }

    if (bottomLeft) {
      drawnCell = this.placeLine(4 * Math.PI / 3, sideLength, drawnCell);
    }

    if (bottomRight) {
      drawnCell = this.placeLine(5 * Math.PI / 3, sideLength, drawnCell);
    }

    // RECALL: 30 degree towards the +Y axis from the +X axis. THIS IS IN GRAPHICS
    // SPACE (not absolute space), so it rotates CLOCKWISE with respect to absolute
    // space
    return new RotateImage(drawnCell, 30);
  }

  // Constructs a boxed square cell with the given background color, wall color,
  // and side length
  private WorldImage boxedImage(int sideLength, Color backgroundColor, Color wallColor) {
    return this.boxedImage2(sideLength, backgroundColor, wallColor, true, true, true, true);
  }

  // Constructs a boxed hexagon cell with the given background color, wall color,
  // and side length
  private WorldImage boxedHexagonImage2(int sideLength, Color backgroundColor, Color wallColor) {
    return this.boxedHexagonImage(sideLength, backgroundColor, wallColor, true, true, true, true,
        true, true);
  }

  // Constructs a boxed square cell with the given background color, wall color,
  // and side length, and
  // whether or not to draw the top, right, left, and bottom walls
  private WorldImage boxedImage2(int sideLength, Color backgroundColor, Color wallColor,
      boolean top, boolean left, boolean right, boolean bottom) {
    WorldImage newBox = new RectangleImage(sideLength, sideLength, OutlineMode.SOLID,
        backgroundColor);
    if (top) {
      WorldImage topWall = new LineImage(new Posn(sideLength, 0), wallColor);
      newBox = new OverlayOffsetImage(topWall, 0, (double) (sideLength) / 2.0, newBox);
    }
    if (left) {
      WorldImage leftWall = new LineImage(new Posn(0, sideLength), wallColor);
      newBox = new OverlayOffsetImage(leftWall, (double) (sideLength) / 2.0, 0, newBox);
    }
    if (right) {
      WorldImage rightWall = new LineImage(new Posn(0, sideLength), wallColor);
      newBox = new OverlayOffsetImage(rightWall, (double) (-sideLength) / 2.0, 0, newBox);
    }
    if (bottom) {
      WorldImage bottomWall = new LineImage(new Posn(sideLength, 0), wallColor);
      newBox = new OverlayOffsetImage(bottomWall, 0.0, (double) (-sideLength) / 2.0, newBox);
    }
    return newBox;
  }

  // A helper method that makes drawing into scenes much faster. Given 9 different
  // colors to draw into each cell of the _threeByThreeHexMaze_, draws the
  // _threeByThreeHexMaze_
  // with the colors in the first row as `col00`, `col01`, etc, the colors in the
  // second row
  // as `col10`, `col11`, `col12`, etc.
  private WorldScene drawExample3x3HexagonMazeWithColors(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(800, 800);

    // First row

    WorldImage hex1 = this.boxedHexagonImage(115, col00, Constants.WALL_COLOR, true, true, true,
        false, true, true);

    drawnManually.placeImageXY(hex1.movePinhole(-hex1.getWidth() / 2, -hex1.getHeight() / 2), 0, 0);

    WorldImage hex2 = this.boxedHexagonImage(115, col01, Constants.WALL_COLOR, true, true, false,
        true, false, true);

    drawnManually.placeImageXY(hex2.movePinhole(-hex2.getWidth() / 2, -hex2.getHeight() / 2), 198,
        0);

    WorldImage hex3 = this.boxedHexagonImage(115, col02, Constants.WALL_COLOR, true, true, true,
        true, true, false);

    drawnManually.placeImageXY(hex3.movePinhole(-hex3.getWidth() / 2, -hex3.getHeight() / 2), 396,
        0);

    // Middle row

    WorldImage hex4 = this.boxedHexagonImage(115, col10, Constants.WALL_COLOR, false, true, true,
        true, true, false);

    drawnManually.placeImageXY(hex4.movePinhole(-hex4.getWidth() / 2, -hex4.getHeight() / 2), 99,
        171);

    WorldImage hex5 = this.boxedHexagonImage(115, col11, Constants.WALL_COLOR, true, true, true,
        false, false, true);

    drawnManually.placeImageXY(hex5.movePinhole(-hex5.getWidth() / 2, -hex5.getHeight() / 2), 297,
        171);

    WorldImage hex6 = this.boxedHexagonImage(115, col12, Constants.WALL_COLOR, true, false, false,
        true, false, true);

    drawnManually.placeImageXY(hex6.movePinhole(-hex6.getWidth() / 2, -hex6.getHeight() / 2), 495,
        171);

    // Bottom row

    WorldImage hex7 = this.boxedHexagonImage(115, col20, Constants.WALL_COLOR, true, true, true,
        false, true, true);

    drawnManually.placeImageXY(hex7.movePinhole(-hex7.getWidth() / 2, -hex7.getHeight() / 2), 0,
        342);

    WorldImage hex8 = this.boxedHexagonImage(115, col21, Constants.WALL_COLOR, false, false, false,
        true, true, true);

    drawnManually.placeImageXY(hex8.movePinhole(-hex8.getWidth() / 2, -hex8.getHeight() / 2), 198,
        342);

    WorldImage hex9 = this.boxedHexagonImage(115, col22, Constants.WALL_COLOR, false, true, true,
        true, true, true);

    drawnManually.placeImageXY(hex9.movePinhole(-hex9.getWidth() / 2, -hex9.getHeight() / 2), 396,
        342);

    return drawnManually;
  }

  // A helper method that makes drawing into scenes much faster. Given 9 different
  // colors to draw into each cell of the _threeByThreeHexMaze_, draws the
  // _threeByThreeHexMaze_
  // with the colors in the first row as `col00`, `col01`, etc, the colors in the
  // second row
  // as `col10`, `col11`, `col12`, etc.
  private WorldScene drawExample3x3HexagonMazeWithColors2(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(90, 90);

    // First row

    WorldImage hex1 = this.boxedHexagonImage(13, col00, Constants.WALL_COLOR, true, true, true,
        false, true, false);

    drawnManually.placeImageXY(hex1.movePinhole(-hex1.getWidth() / 2, -hex1.getHeight() / 2), 0, 0);

    WorldImage hex2 = this.boxedHexagonImage(13, col01, Constants.WALL_COLOR, true, true, false,
        true, true, false);

    drawnManually.placeImageXY(hex2.movePinhole(-hex2.getWidth() / 2, -hex2.getHeight() / 2), 22,
        0);

    WorldImage hex3 = this.boxedHexagonImage(13, col02, Constants.WALL_COLOR, true, true, true,
        true, false, true);

    drawnManually.placeImageXY(hex3.movePinhole(-hex3.getWidth() / 2, -hex3.getHeight() / 2), 44,
        0);

    // Middle row

    WorldImage hex4 = this.boxedHexagonImage(13, col10, Constants.WALL_COLOR, true, false, true,
        true, false, true);

    drawnManually.placeImageXY(hex4.movePinhole(-hex4.getWidth() / 2, -hex4.getHeight() / 2), 11,
        18);

    WorldImage hex5 = this.boxedHexagonImage(13, col11, Constants.WALL_COLOR, false, false, true,
        false, true, true);

    drawnManually.placeImageXY(hex5.movePinhole(-hex5.getWidth() / 2, -hex5.getHeight() / 2), 33,
        18);

    WorldImage hex6 = this.boxedHexagonImage(13, col12, Constants.WALL_COLOR, true, true, false,
        true, true, true);

    drawnManually.placeImageXY(hex6.movePinhole(-hex6.getWidth() / 2, -hex6.getHeight() / 2), 55,
        18);

    // Bottom row

    WorldImage hex7 = this.boxedHexagonImage(13, col20, Constants.WALL_COLOR, false, true, true,
        false, true, true);

    drawnManually.placeImageXY(hex7.movePinhole(-hex7.getWidth() / 2, -hex7.getHeight() / 2), 0,
        36);

    WorldImage hex8 = this.boxedHexagonImage(13, col21, Constants.WALL_COLOR, true, true, false,
        false, true, true);

    drawnManually.placeImageXY(hex8.movePinhole(-hex8.getWidth() / 2, -hex8.getHeight() / 2), 22,
        36);

    WorldImage hex9 = this.boxedHexagonImage(13, col22, Constants.WALL_COLOR, true, true, false,
        true, true, true);

    drawnManually.placeImageXY(hex9.movePinhole(-hex9.getWidth() / 2, -hex9.getHeight() / 2), 44,
        36);

    return drawnManually;
  }

  // A helper method that makes drawing into scenes much faster. Given 9 different
  // colors to draw into each cell of the _threeByThreeMaze_, draws the
  // _threeByThreeMaze_
  // with the colors in the first row as `col00`, `col01`, etc, the colors in the
  // second row
  // as `col10`, `col11`, `col12`, etc.
  private WorldScene drawExample3x3MazeWithColors(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(90, 90);

    // First row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col00, Constants.WALL_COLOR, 
            true, true, true, false).movePinhole(-15, -15), 0, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col01, Constants.WALL_COLOR, 
            true, true, false, true).movePinhole(-15, -15), 30, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col02, Constants.WALL_COLOR, 
            true, false, true, false).movePinhole(-15, -15), 60, 0);

    // Middle row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col10, Constants.WALL_COLOR, false, true, false, false)
            .movePinhole(-15, -15), 0, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col11, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col12, Constants.WALL_COLOR, false, false, true, true)
            .movePinhole(-15, -15), 60, 30);

    // Bottom row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col20, Constants.WALL_COLOR, false, true, false, true)
            .movePinhole(-15, -15), 0, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col21, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col22, Constants.WALL_COLOR, 
            true, false, true, true).movePinhole(-15, -15), 60, 60);

    return drawnManually;
  }

  // makes a maze in the same way as drawExample3x3MazeWithColors, just making the
  // first maze that is created after making a new unbiased maze on the prior
  private WorldScene drawExample3x3MazeWithColors2(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(90, 90);

    // First row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col00, Constants.WALL_COLOR, true, true, false, false)
            .movePinhole(-15, -15), 0, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col01, Constants.WALL_COLOR, 
            true, false, true, true).movePinhole(-15, -15), 30, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col02, Constants.WALL_COLOR, 
            true, true, true, false).movePinhole(-15, -15), 60, 0);

    // Middle row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col10, Constants.WALL_COLOR, false, true, false, true)
            .movePinhole(-15, -15), 0, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col11, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col12, Constants.WALL_COLOR, false, false, true, false)
            .movePinhole(-15, -15), 60, 30);

    // Bottom row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col20, Constants.WALL_COLOR, 
            true, true, false, true).movePinhole(-15, -15), 0, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col21, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col22, Constants.WALL_COLOR, false, false, true, true)
            .movePinhole(-15, -15), 60, 60);

    return drawnManually;
  }

  // makes a maze in the same way as drawExample3x3MazeWithColors2, just making
  // the
  // first maze that is created after making a new rowed maze on the prior
  private WorldScene drawExample3x3MazeWithColorsRow(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(90, 90);

    // First row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col00, Constants.WALL_COLOR, 
            true, true, false, true).movePinhole(-15, -15), 0, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col01, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col02, Constants.WALL_COLOR, true, false, true, false)
            .movePinhole(-15, -15), 60, 0);

    // Middle row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col10, Constants.WALL_COLOR, true, true, false, false)
            .movePinhole(-15, -15), 0, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col11, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col12, Constants.WALL_COLOR, false, false, true, true)
            .movePinhole(-15, -15), 60, 30);

    // Bottom row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col20, Constants.WALL_COLOR, false, true, false, true)
            .movePinhole(-15, -15), 0, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col21, Constants.WALL_COLOR, true, false, false, true)
            .movePinhole(-15, -15), 30, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col22, Constants.WALL_COLOR, 
            true, false, true, true).movePinhole(-15, -15), 60, 60);

    return drawnManually;
  }

  // makes a maze in the same way as drawExample3x3MazeWithColorsRow, just making
  // the
  // first maze that is created after making a new columned maze on the prior
  private WorldScene drawExample3x3MazeWithColorsCol(Color col00, Color col01, Color col02,
      Color col10, Color col11, Color col12, Color col20, Color col21, Color col22) {
    WorldScene drawnManually = new WorldScene(90, 90);

    // First row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col00, Constants.WALL_COLOR, true, true, false, false)
            .movePinhole(-15, -15), 0, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col01, Constants.WALL_COLOR, true, false, false, false)
            .movePinhole(-15, -15), 30, 0);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col02, Constants.WALL_COLOR, true, false, true, false)
            .movePinhole(-15, -15), 60, 0);

    // Middle row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col10, Constants.WALL_COLOR, false, true, true, false)
            .movePinhole(-15, -15), 0, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col11, Constants.WALL_COLOR, false, true, true, false)
            .movePinhole(-15, -15), 30, 30);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col12, Constants.WALL_COLOR, false, true, true, false)
            .movePinhole(-15, -15), 60, 30);

    // Bottom row
    drawnManually
        .placeImageXY(this.boxedImage2(30, col20, Constants.WALL_COLOR, 
            false, true, true, true).movePinhole(-15, -15), 0, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col21, Constants.WALL_COLOR, 
            false, true, true, true).movePinhole(-15, -15), 30, 60);
    drawnManually
        .placeImageXY(this.boxedImage2(30, col22, Constants.WALL_COLOR, 
            false, true, true, true).movePinhole(-15, -15), 60, 60);

    return drawnManually;
  }

  void initTestData() {

    // The value `-1155484576` is the first random integer generated
    // by `new Random(0)`. The value of `0` given as the `seedGenerationSeed`
    // of the `MazeWorld` constructor would produce as its first Maze the
    // Maze with this seed
    int firstIntForRandom0 = -1155484576;
    this.oneByOneMaze = new Maze(1, 1, 10);
    this.twoByTwoMaze = new Maze(2, 2, firstIntForRandom0);
    this.threeByThreeMaze = new Maze(3, 3, firstIntForRandom0);
    this.threeByThreeHexMaze = new Maze(3, 3, firstIntForRandom0, true);

    // Knock down all of the walls for drawing
    this.threeByThreeMaze.breakAllWalls();
    this.twoByTwoMaze.breakAllWalls();
    this.threeByThreeMaze.breakAllWalls();
    this.threeByThreeHexMaze.breakAllWalls();

    this.boxCellImageUnseen100x100 = this.boxedImage(100, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageSeen100x100 = this.boxedImage(100, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageStart100x100 = this.boxedImage(100, Constants.CELL_START_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageFinish100x100 = this.boxedImage(100, Constants.CELL_FINISH_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImagePath100x100 = this.boxedImage(100, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR);

    this.boxCellImageUnseenInterpolated100x100 = this.boxedImage(100,
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImageSeenInterpolated100x100 = this.boxCellImageSeen100x100;
    this.boxCellImageStartInterpolated100x100 = this.boxedImage(100,
        new ColorUtils().interpolate(0.5, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImageFinishInterpolated100x100 = this.boxedImage(100,
        new ColorUtils().interpolate(0.5, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImagePathInterpolated100x100 = this.boxCellImagePath100x100;

    this.hexCellImageUnseenInterpolated100x100 = this.boxedHexagonImage2(100,
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.hexCellImageSeenInterpolated100x100 = this.hexCellImageSeen100x100;
    this.hexCellImageStartInterpolated100x100 = this.boxedHexagonImage2(100,
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.hexCellImageFinishInterpolated100x100 = this.boxedHexagonImage2(100,
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.hexCellImagePathInterpolated100x100 = this.hexCellImagePath100x100;

    this.boxCellImageUnseen150x150 = this.boxedImage(150, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageSeen150x150 = this.boxedImage(150, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageStart150x150 = this.boxedImage(150, Constants.CELL_START_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImageFinish150x150 = this.boxedImage(150, Constants.CELL_FINISH_COLOR,
        Constants.WALL_COLOR);
    this.boxCellImagePath150x150 = this.boxedImage(150, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR);

    this.boxCellImageUnseenInterpolated150x150 = this.boxedImage(150,
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImageSeenInterpolated150x150 = this.boxCellImageSeen150x150;
    this.boxCellImageStartInterpolated150x150 = this.boxedImage(150,
        new ColorUtils().interpolate(0.5, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImageFinishInterpolated150x150 = this.boxedImage(150,
        new ColorUtils().interpolate(0.5, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR),
        Constants.WALL_COLOR);
    this.boxCellImagePathInterpolated150x150 = this.boxCellImagePath150x150;

    // Images of hexagonal cells of size 100x100 with all of its walls
    // in the various states
    this.hexCellImageUnseen100x100 = this.boxedHexagonImage2(100, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR);
    this.hexCellImageSeen100x100 = this.boxedHexagonImage2(100, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR);
    this.hexCellImageStart100x100 = this.boxedHexagonImage2(100, Constants.CELL_START_COLOR,
        Constants.WALL_COLOR);
    this.hexCellImageFinish100x100 = this.boxedHexagonImage2(100, Constants.CELL_FINISH_COLOR,
        Constants.WALL_COLOR);
    this.hexCellImagePath100x100 = this.boxedHexagonImage2(100, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR);
  }

  void testBigBang(Tester t) {
    this.initTestData();

    MazeWorld w = new MazeWorld(200, 200, Constants.WINDOW_WIDTH, Constants.WINDOW_HEIGHT, 0,
        false);
    w.bigBang(Constants.WINDOW_WIDTH, Constants.WINDOW_HEIGHT, .01);
  }

  void testConnectAtRightAndBottom(Tester t) {
    this.initTestData();

    Cell sampleCellLeft = new Cell();
    Cell sampleCellRight = new Cell();
    Cell sampleCellTop = new Cell();
    Cell sampleCellBottom = new Cell();
    Cell sampleCellCenter = new Cell();

    // begin by checking that we are getting the appropriate edges out of
    // our connection methods
    t.checkExpect(sampleCellLeft.connectAtRight(sampleCellCenter, 10),
        new CellEdge(sampleCellLeft, sampleCellCenter, 10));
    t.checkExpect(sampleCellCenter.connectAtRight(sampleCellRight, 20),
        new CellEdge(sampleCellCenter, sampleCellRight, 20));

    t.checkExpect(sampleCellTop.connectAtBottom(sampleCellCenter, 30),
        new CellEdge(sampleCellTop, sampleCellCenter, 30));
    t.checkExpect(sampleCellCenter.connectAtBottom(sampleCellBottom, 40),
        new CellEdge(sampleCellCenter, sampleCellBottom, 40));

    // we can then confirm that our connections are in place by testing the cells
    // that are returned when using the allAdjacent
    t.checkExpect(sampleCellCenter.allAdjacentCells(), new ArrayList<ICell>(
        List.of(sampleCellTop, sampleCellLeft, sampleCellRight, sampleCellBottom)));

    // and we'll do another
    Cell sampleCellLeft2 = new Cell();
    Cell sampleCellRight2 = new Cell();
    Cell sampleCellTop2 = new Cell();
    Cell sampleCellBottom2 = new Cell();
    Cell sampleCellCenter2 = new Cell();

    t.checkExpect(sampleCellLeft2.connectAtRight(sampleCellCenter2, 0),
        new CellEdge(sampleCellLeft2, sampleCellCenter2, 0));
    t.checkExpect(sampleCellCenter2.connectAtRight(sampleCellRight2, 0),
        new CellEdge(sampleCellCenter2, sampleCellRight2, 0));

    t.checkExpect(sampleCellTop2.connectAtBottom(sampleCellCenter2, 0),
        new CellEdge(sampleCellTop2, sampleCellCenter2, 0));
    t.checkExpect(sampleCellCenter2.connectAtBottom(sampleCellBottom2, 0),
        new CellEdge(sampleCellCenter2, sampleCellBottom2, 0));

    t.checkExpect(sampleCellCenter2.allAdjacentCells(), new ArrayList<ICell>(
        List.of(sampleCellTop2, sampleCellLeft2, sampleCellRight2, sampleCellBottom2)));

  }

  void testMazeConstructor(Tester t) {
    this.initTestData();

    // Ensure that constructing a maze with a
    // negative number of cells per row or per column (or 0 per row/column)
    // doesn't make sense
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per row. Given -10"),
        "Maze", -10, 11, 0);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per column. Given -4"),
        "Maze", 5, -4, -190);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per row. Given 0"), "Maze",
        0, 5, 0);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per column. Given 0"),
        "Maze", 5, 0, -4);

    // It doesn't make sense to have both horizontal and vertical bias
    t.checkConstructorException(
        new IllegalArgumentException(
            "A maze cannot be both horizontally and vertically biased at the same time"),
        "Maze", 5, 6, 100, true, true, false);
  }

  void testMazeWorldConstructor(Tester t) {
    this.initTestData();

    // Ensure that constructing a maze with a
    // negative number of cells per row or per column (or 0 per row/column)
    // doesn't make sense
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per row. Given -10"),
        "MazeWorld", -10, 11, 100, 100, 0, false);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per column. Given -4"),
        "MazeWorld", 5, -4, 100, 100, -190, false);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per row. Given 0"),
        "MazeWorld", 0, 5, 100, 100, 0, false);
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have at least one cell per column. Given 0"),
        "MazeWorld", 5, 0, 100, 100, -4, false);

    // It doesn't make sense to have a negative window width or height
    t.checkConstructorException(
        new IllegalArgumentException("A maze must have a positive window width"), "MazeWorld", 5, 5,
        -100, 100, -4, false);

    t.checkConstructorException(
        new IllegalArgumentException("A maze must have a positive window height"), "MazeWorld", 5,
        5, 100, -100, -4, false);
  }

  void testWallMethods(Tester t) {

    // tests for breakAllWalls, breakOneWall, and doneBreakingWalls

    Maze smallMaze = new Maze(3, 3, 7, false, false, false);
    // this maze has 12 walls to start, and is left with 4
    // after all walls are broken

    t.checkExpect(smallMaze.doneBreakingWalls(), false);
    smallMaze.breakOneWall();
    smallMaze.breakOneWall();
    smallMaze.breakOneWall();
    smallMaze.breakOneWall();

    // still has 4 more walls to break
    t.checkExpect(smallMaze.doneBreakingWalls(), false);
    smallMaze.breakOneWall();
    smallMaze.breakOneWall();

    // still 2 more to break
    t.checkExpect(smallMaze.doneBreakingWalls(), false);

    smallMaze.breakOneWall();
    smallMaze.breakOneWall();
    t.checkExpect(smallMaze.doneBreakingWalls(), true);

    // we don't really care about this mazes wall count, as we are
    // going to break them all at once
    Maze smallMaze2 = new Maze(3, 3, 6, false, false, false);

    smallMaze2.breakAllWalls();
    t.checkExpect(smallMaze2.doneBreakingWalls(), true);

    // here we will break some walls one at a time, then clear
    // the rest of them all at once
    Maze smallMaze3 = new Maze(3, 3, 6, false, false, false);

    smallMaze3.breakOneWall();
    smallMaze3.breakOneWall();
    smallMaze3.breakOneWall();
    smallMaze3.breakOneWall();

    t.checkExpect(smallMaze3.doneBreakingWalls(), false);

    // now we break the rest
    smallMaze3.breakAllWalls();
    t.checkExpect(smallMaze3.doneBreakingWalls(), true);
  }

  void testDrawInScene(Tester t) {
    this.initTestData();

    // 3x3 hexagon maze
    WorldScene drawnByMaze = new WorldScene(800, 800);
    WorldScene drawnManually = new WorldScene(800, 800);

    // Drawn using the Maze
    this.threeByThreeHexMaze.drawInScene(drawnByMaze, 800, 800);

    drawnManually = this.drawExample3x3HexagonMazeWithColors(Constants.CELL_START_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    // Should match
    t.checkExpect(drawnByMaze, drawnManually);

    // 3x3 regular maze
    drawnByMaze = new WorldScene(90, 90);
    drawnManually = new WorldScene(90, 90);

    // Drawn using the Maze
    this.threeByThreeMaze.drawInScene(drawnByMaze, 90, 90);

    // Manually drawn
    // The maze drawn out looks like
    // __ __ __
    // | |__ |
    // | _____|
    // |__ __ __|
    //

    drawnManually = this.drawExample3x3MazeWithColors(Constants.CELL_START_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    // Should match
    t.checkExpect(drawnByMaze, drawnManually);

    // 2x2 maze
    drawnByMaze = new WorldScene(60, 60);
    drawnManually = new WorldScene(60, 60);

    this.twoByTwoMaze.drawInScene(drawnByMaze, 60, 60);

    // Manually drawn
    // The maze drawn out looks like
    // __ __
    // | __|
    // |__ __|
    //

    // First row
    drawnManually.placeImageXY(this
        .boxedImage2(30, Constants.CELL_START_COLOR, Constants.WALL_COLOR, true, true, false, false)
        .movePinhole(-15, -15), 0, 0);
    drawnManually.placeImageXY(this
        .boxedImage2(30, Constants.CELL_UNSEEN_COLOR, Constants.WALL_COLOR, true, false, true, true)
        .movePinhole(-15, -15), 30, 0);

    // Second row
    drawnManually.placeImageXY(this.boxedImage2(30, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, false, true).movePinhole(-15, -15), 0, 30);
    drawnManually.placeImageXY(this
        .boxedImage2(30, Constants.CELL_FINISH_COLOR, Constants.WALL_COLOR, 
            true, false, true, true).movePinhole(-15, -15), 30, 30);

    t.checkExpect(drawnByMaze, drawnManually);

    // 1x1 maze
    drawnByMaze = new WorldScene(30, 30);
    drawnManually = new WorldScene(30, 30);

    this.oneByOneMaze.drawInScene(drawnByMaze, 30, 30);

    drawnManually = new WorldScene(30, 30);

    drawnManually.placeImageXY(this
        .boxedImage2(30, Constants.CELL_FINISH_COLOR, Constants.WALL_COLOR, 
            true, true, true, true) .movePinhole(-15, -15), 0, 0);

    t.checkExpect(drawnByMaze, drawnManually);
  }

  // Produces the color of a cell a distance _distance_ away from the start (or
  // finish)
  // given a maximum distance from the start or finish
  private Color cellColorInterpolated(int distance, int max) {
    return new ColorUtils().interpolate((double) (distance) / (double) (max),
        Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR);
  }

  void testDrawInSceneInterpolated(Tester t) {
    this.initTestData();

    // 3x3 maze
    WorldScene drawnByMaze = new WorldScene(90, 90);
    WorldScene drawnManually = new WorldScene(90, 90);

    // Drawn using the Maze
    this.threeByThreeMaze.drawInSceneInterpolated(drawnByMaze, 90, 90, true);

    // Manually drawn
    // The maze drawn out looks like
    // __ __ __
    // | |_*_ |
    // | _____|
    // |__ _ _*_|
    //
    // Max distance from start: 5 (shown as *)

    drawnManually = this.drawExample3x3MazeWithColors(this.cellColorInterpolated(0, 5),
        this.cellColorInterpolated(5, 5), this.cellColorInterpolated(4, 5),
        this.cellColorInterpolated(1, 5), this.cellColorInterpolated(2, 5),
        this.cellColorInterpolated(3, 5), this.cellColorInterpolated(2, 5),
        this.cellColorInterpolated(3, 5), this.cellColorInterpolated(4, 5));

    // Should match
    t.checkExpect(drawnByMaze, drawnManually);

    // 3x3 maze, interpolation from the final cell
    drawnByMaze = new WorldScene(90, 90);
    drawnManually = new WorldScene(90, 90);

    // Drawn using the Maze
    this.threeByThreeMaze.drawInSceneInterpolated(drawnByMaze, 90, 90, false);

    // Manually drawn
    // The maze drawn out looks like
    // __ __ __
    // | |_*_ |
    // | _____|
    // |__ __ __|
    //
    // Max distance from end: 7 (shown in *)

    drawnManually = this.drawExample3x3MazeWithColors(this.cellColorInterpolated(4, 7),
        this.cellColorInterpolated(7, 7), this.cellColorInterpolated(6, 7),
        this.cellColorInterpolated(3, 7), this.cellColorInterpolated(4, 7),
        this.cellColorInterpolated(5, 7), this.cellColorInterpolated(2, 7),
        this.cellColorInterpolated(1, 7), this.cellColorInterpolated(0, 7));

    // Should match
    t.checkExpect(drawnByMaze, drawnManually);

    // 2x2 maze
    drawnByMaze = new WorldScene(60, 60);
    drawnManually = new WorldScene(60, 60);

    this.twoByTwoMaze.drawInSceneInterpolated(drawnByMaze, 60, 60, true);

    // Manually drawn
    // The maze drawn out looks like
    // __ __
    // | __|
    // |__ __|
    //
    // Max distance: 2

    // First row
    drawnManually.placeImageXY(this.boxedImage2(30, this.cellColorInterpolated(0, 2),
        Constants.WALL_COLOR, true, true, false, false).movePinhole(-15, -15), 0, 0);
    drawnManually.placeImageXY(this.boxedImage2(30, this.cellColorInterpolated(1, 2),
        Constants.WALL_COLOR, true, false, true, true).movePinhole(-15, -15), 30, 0);

    // Second row
    drawnManually.placeImageXY(this.boxedImage2(30, this.cellColorInterpolated(1, 2),
        Constants.WALL_COLOR, false, true, false, true).movePinhole(-15, -15), 0, 30);
    drawnManually.placeImageXY(this.boxedImage2(30, this.cellColorInterpolated(2, 2),
        Constants.WALL_COLOR, true, false, true, true).movePinhole(-15, -15), 30, 30);

    t.checkExpect(drawnByMaze, drawnManually);

    // 1x1 maze
    drawnByMaze = new WorldScene(30, 30);
    drawnManually = new WorldScene(30, 30);

    this.oneByOneMaze.drawInSceneInterpolated(drawnByMaze, 30, 30, true);

    drawnManually = new WorldScene(30, 30);

    drawnManually.placeImageXY(this.boxedImage2(30, this.cellColorInterpolated(0, 0),
        Constants.WALL_COLOR, true, true, true, true).movePinhole(-15, -15), 0, 0);

    t.checkExpect(drawnByMaze, drawnManually);
  }

  // The way we test the `new___Snapshot()` methods on
  // Maze (shown below) is as follows:
  //
  // 1. Create the snapshot we are trying to test
  // 2. Have it change the state of the cells it has
  // run into using `makeAllSeen()` and then have the
  // maze draw into an empty scene. This ensures
  // that we'll be able to see the progress of the snapshot
  // as we move it through the Maze (see the next step)
  // 3. Call `searchOneStep()` on the snapshot. Again, call
  // `makeAllSeen()` on the snapshot and again have
  // the maze draw into another empty scene.
  // 4. Repeat until the search has finished
  //
  // This method allows us to incrementally inspect the states
  // of the cells in the maze as well as the cells that the snapshot has
  // processed without needing access to the cells themselves. We are effectively
  // looking
  // to test the behavior of the snapshots through the images it produces in
  // scenes

  void testNewDepthFirstSnapshot(Tester t) {
    this.initTestData();

    WorldScene drawnScene = new WorldScene(90, 90);
    ASearchSnapshot depthSnapshot = this.threeByThreeMaze.newDepthFirstSnapshot();

    // INITIAL

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Draw out what we think the scene should be at this point
    WorldScene expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_START_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 1:

    // Now, move the snapshot through the scene
    boolean doneNow = depthSnapshot.searchOneStep();
    depthSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 2:

    // Move another step
    doneNow = depthSnapshot.searchOneStep();
    depthSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 3:

    // Move the snapshot through the scene AGAIN
    doneNow = depthSnapshot.searchOneStep();
    depthSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene yet again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 4:

    // Move the snapshot through the scene AGAIN
    doneNow = depthSnapshot.searchOneStep();
    depthSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene yet again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_HEAD_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 5:

    // Move the snapshot through the scene AGAIN
    doneNow = depthSnapshot.searchOneStep();
    depthSnapshot.makeAllSeen();

    // We have hit the end
    t.checkExpect(doneNow, true);

    // Draw the scene yet again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_HEAD_COLOR);

    t.checkExpect(drawnScene, expectedScene);

  }

  void testNewBreadthFirstSnapshot(Tester t) {
    this.initTestData();

    WorldScene drawnScene = new WorldScene(90, 90);
    ASearchSnapshot breadthFirstSnapshot = this.threeByThreeMaze.newBreadthFirstSnapshot();

    // INITIAL

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Draw out what we think the scene should be at this point
    WorldScene expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_START_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 1:

    // Now, move the snapshot through the scene
    boolean doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 2:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 3:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_HEAD_COLOR, Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 4:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 5:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 6:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_HEAD_COLOR, Constants.CELL_FUTUREHEAD_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 7:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    t.checkExpect(doneNow, false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // STEP 8:

    // Move another step
    doneNow = breadthFirstSnapshot.searchOneStep();
    breadthFirstSnapshot.makeAllSeen();

    // We are now done
    t.checkExpect(doneNow, true);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // Using depth first, we move as far as possible
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_FUTUREHEAD_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_HEAD_COLOR);

    t.checkExpect(drawnScene, expectedScene);
  }

  void testNewExhaustedSnapshotFromStart(Tester t) {
    this.initTestData();

    WorldScene drawnScene = new WorldScene(90, 90);
    ASearchSnapshot exhaustedSnapshot = this.threeByThreeMaze.newExhaustedSnapshotFromStart();

    // To test an exhausted snapshot
    exhaustedSnapshot.makeAllSeen();

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    WorldScene expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_HEAD_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // To ensure that the depth snapshot correctly
    // captured the distances of all cells in the maze,
    // draw the scene using interpolated colors. Thi

    // Also, check that the farthest distance is as we expect
    t.checkExpect(exhaustedSnapshot.distanceToFarthest(), 5);

    // Unsee everything to test color interpolation
    exhaustedSnapshot.makeAllUnseen();

    t.checkExpect(exhaustedSnapshot.hasMoreToSee(), false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInSceneInterpolated(drawnScene, 90, 90, true);

    // In order to test that drawing an exhausted search snapshot works
    // as we expect it to, we test that the interpolated color image
    // is produced correctly, since it uses such a snapshot itself

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(this.cellColorInterpolated(0, 5),
        this.cellColorInterpolated(5, 5), this.cellColorInterpolated(4, 5),
        this.cellColorInterpolated(1, 5), this.cellColorInterpolated(2, 5),
        this.cellColorInterpolated(3, 5), this.cellColorInterpolated(2, 5),
        this.cellColorInterpolated(3, 5), this.cellColorInterpolated(4, 5));

    t.checkExpect(drawnScene, expectedScene);
  }

  void testNewExhaustedSnapshotFromFinish(Tester t) {
    this.initTestData();

    WorldScene drawnScene = new WorldScene(90, 90);
    ASearchSnapshot exhaustedSnapshot = this.threeByThreeMaze.newExhaustedSnapshotFromFinish();

    // To test an exhausted snapshot
    exhaustedSnapshot.makeAllSeen();

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    WorldScene expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // To ensure that the depth snapshot correctly
    // captured the distances of all cells in the maze,
    // draw the scene using interpolated colors. Thi

    // Also, check that the farthest distance is as we expect
    t.checkExpect(exhaustedSnapshot.distanceToFarthest(), 7);

    // Unsee everything to test color interpolation
    exhaustedSnapshot.makeAllUnseen();

    t.checkExpect(exhaustedSnapshot.hasMoreToSee(), false);

    // Draw the scene again
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInSceneInterpolated(drawnScene, 90, 90, false);

    // In order to test that drawing an exhausted search snapshot works
    // as we expect it to, we test that the interpolated color image
    // is produced correctly, since it uses such a snapshot itself

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(this.cellColorInterpolated(4, 7),
        this.cellColorInterpolated(7, 7), this.cellColorInterpolated(6, 7),
        this.cellColorInterpolated(3, 7), this.cellColorInterpolated(4, 7),
        this.cellColorInterpolated(5, 7), this.cellColorInterpolated(2, 7),
        this.cellColorInterpolated(1, 7), this.cellColorInterpolated(0, 7));

    t.checkExpect(drawnScene, expectedScene);
  }

  void testNewManualSnapshot(Tester t) {
    this.initTestData();

    WorldScene drawnScene = new WorldScene(90, 90);
    ASearchSnapshot manualSnapshot = this.threeByThreeMaze.newManualSnapshot();

    // To test that initially the player
    // starts at the beginning
    manualSnapshot.makeAllSeen();

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    WorldScene expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move down
    // We know this will move down in our maze
    manualSnapshot.maybeMoveDown();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move down AGAIN
    // We know this will move down in our maze
    manualSnapshot.maybeMoveDown();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move down AGAIN. Nothing will happen in this case
    manualSnapshot.maybeMoveDown();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move UP. The path
    // the player has travelled will be left behind
    manualSnapshot.maybeMoveUp();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move RIGHT. The path
    // the player has travelled will be left behind
    manualSnapshot.maybeMoveRight();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);

    // Have the player move LEFT. The path
    // the player has travelled will be left behind
    manualSnapshot.maybeMoveLeft();
    manualSnapshot.makeAllSeen();

    // Create a fresh scene to draw into
    drawnScene = new WorldScene(90, 90);

    // Draw the maze into an empty scene
    this.threeByThreeMaze.drawInScene(drawnScene, 90, 90);

    // At this point, we expected that the exhausted snapshot
    expectedScene = this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
        Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR,
        Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR);

    t.checkExpect(drawnScene, expectedScene);
  }

  void testEdgeConstructor(Tester t) {
    this.initTestData();
    ICell icell1 = new Cell();
    ICell icell2 = new DummyCell();
    ICell icell3 = new Cell();
    ICell icell4 = new DummyCell();
    ICell icell5 = new Cell();
    CellEdge edge1 = new CellEdge(icell1, icell2, 10);
    CellEdge edge2 = new CellEdge(icell3, icell4, 4);
    CellEdge edge3 = new CellEdge(icell4, icell5, 100);
    CellEdge edge4 = new CellEdge(icell2, icell2, 19);
    // Check that the edge that was returned has the correct references and has the
    // correct weight
    t.checkExpect(edge1.first, icell1);
    t.checkExpect(edge1.second, icell2);
    t.checkExpect(edge1.weight, 10);
    t.checkExpect(edge2.first, icell3);
    t.checkExpect(edge2.second, icell4);
    t.checkExpect(edge2.weight, 4);
    t.checkExpect(edge3.first, icell4);
    t.checkExpect(edge3.second, icell5);
    t.checkExpect(edge3.weight, 100);
    t.checkExpect(edge4.first, icell2);
    t.checkExpect(edge4.second, icell2);
    t.checkExpect(edge4.weight, 19);
    // Attempting to make an edge with a negative weight is
    // disallowed in our program (Kruskal's algorithm only likes positive
    // weights)
    t.checkConstructorException(
        new IllegalArgumentException(
            "Cannot represent a connection between two ICells with a negative weight"),
        "CellEdge", icell1, icell2, -10);
    t.checkConstructorException(
        new IllegalArgumentException(
            "Cannot represent a connection between two ICells with a negative weight"),
        "CellEdge", icell4, icell5, -1);
  }

  void testAllPathways(Tester t) {

    // for getting all pathways, we'll make a 3x3 grid that looks like this

    /*   _____
     *  |     |
     *  | |_|_|
     *  |     |
     *   -----
     *
     * each of the walls that remain will be for edges with weight 100,
     * all open edges will have weight 4
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // our expected edges in our minimum spanning tree, based on the weights above
    ArrayList<CellEdge> expectedOpenEdges = new ArrayList<CellEdge>(
        List.of(edge12, edge23, edge78, edge89, edge14, edge25, edge36, edge47));

    t.checkExpect(openEdges, expectedOpenEdges);
  }

  void testUnionFindContains(Tester t) {
    this.initTestData();

    // Construct a new UnionFind structure of strings
    UnionFind<String> stringUnionFind = new UnionFind<String>(
        new ArrayList<String>(List.of("A", "B", "C", "D", "E")));

    // Ensure that the union find contains the given values
    t.checkExpect(stringUnionFind.contains("A"), true);
    t.checkExpect(stringUnionFind.contains("B"), true);
    t.checkExpect(stringUnionFind.contains("C"), true);
    t.checkExpect(stringUnionFind.contains("D"), true);
    t.checkExpect(stringUnionFind.contains("E"), true);

    // It does not contain these values though
    t.checkExpect(stringUnionFind.contains("This value"), false);
    t.checkExpect(stringUnionFind.contains("hello"), false);
  }

  void testUnionFindBehavior(Tester t) {
    this.initTestData();

    // Construct a new UnionFind structure of strings.
    UnionFind<String> stringUnionFind = new UnionFind<String>(
        new ArrayList<String>(List.of("A", "B", "C", "D", "E")));

    // Check that none of the nodes are connected yet
    t.checkExpect(stringUnionFind.findRepresentative("A"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("B"), "B");
    t.checkExpect(stringUnionFind.findRepresentative("C"), "C");
    t.checkExpect(stringUnionFind.findRepresentative("D"), "D");
    t.checkExpect(stringUnionFind.findRepresentative("E"), "E");

    // Union together the trees that these elements are apart of
    stringUnionFind.union("A", "B");

    // Check that the representative of "A" is now "B", and that "A" is still "A"
    t.checkExpect(stringUnionFind.findRepresentative("B"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("A"), "A");

    // Check that the other nodes remained in-tact (that is, that
    // we still have the forest as we expected)
    t.checkExpect(stringUnionFind.findRepresentative("C"), "C");
    t.checkExpect(stringUnionFind.findRepresentative("D"), "D");
    t.checkExpect(stringUnionFind.findRepresentative("E"), "E");

    // Union "C" with the "A" <- "B" tree to make a new tree
    // "A"
    // "B" and "C" point to "A"
    stringUnionFind.union("C", "A");
    t.checkExpect(stringUnionFind.findRepresentative("A"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("B"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("C"), "A");

    // Check that the other nodes remained in-tact (that is, that
    // we still have the forest as we expected)
    t.checkExpect(stringUnionFind.findRepresentative("D"), "D");
    t.checkExpect(stringUnionFind.findRepresentative("E"), "E");

    // Now, combine "D" and "E" into a tree and then add that tree to the
    // "A" structure. We use "E" instead of "D" here, but that shouldn't matter
    // because the implementation is supposed to use the representatives in making
    // unions
    stringUnionFind.union("D", "E");
    stringUnionFind.union("A", "E");

    // Now, check that the forest
    t.checkExpect(stringUnionFind.findRepresentative("A"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("B"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("C"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("D"), "A");
    t.checkExpect(stringUnionFind.findRepresentative("E"), "A");

    // If we now attempt to union an element not in the UnionFind, we expect
    // an exception
    t.checkException(new IllegalArgumentException("Could not perform a union for the"
        + " first provided item as it is not in the set of items visible to this "
        + "UnionFind instance"), stringUnionFind, "union", "T", "A");
    t.checkException(new IllegalArgumentException("Could not perform a union for the"
        + " second provided item as it is not in the set of items visible to this "
        + "UnionFind instance"), stringUnionFind, "union", "A", "H");

    // Test that attempting to find the representative
    // for a node that is not in the UnionFind structure is an exception
    t.checkException(new IllegalArgumentException("Could not locate a representative for the"
        + " provided item as it is not in the set of items visible to this "
        + "UnionFind instance"), stringUnionFind, "findRepresentative", "H");
    t.checkException(new IllegalArgumentException("Could not locate a representative for the"
        + " provided item as it is not in the set of items visible to this "
        + "UnionFind instance"), stringUnionFind, "findRepresentative", "Y");

    // Let's try another union find. This time, we define locals
    // to keep their intentional identities (we didn't need to do so for
    // the strings above because of the way the compiler treats string literals)
    Posn posnA = new Posn(10, 10);
    Posn posnB = new Posn(-10, 10);
    Posn posnC = new Posn(10, -10);
    Posn posnD = new Posn(-10, -10);

    // Construct a new UnionFind structure of Posn objects
    UnionFind<Posn> posnUnionFind = new UnionFind<Posn>(
        new ArrayList<Posn>(List.of(posnA, posnB, posnC, posnD)));

    // Check that none of the nodes are connected yet
    t.checkExpect(posnUnionFind.findRepresentative(posnA), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnB), posnB);
    t.checkExpect(posnUnionFind.findRepresentative(posnC), posnC);
    t.checkExpect(posnUnionFind.findRepresentative(posnD), posnD);

    // Check that Posns that are `.equals()` to each other are not
    // considered to be contained in the UnionFind since it uses intentional
    // equality
    Posn posAEquals = new Posn(10, 10);
    t.checkExpect(posAEquals.equals(posnA), true);
    t.checkExpect(posnUnionFind.contains(posAEquals), false);

    // Connect posnA with posnB
    posnUnionFind.union(posnA, posnB);

    // Check the new representatives
    t.checkExpect(posnUnionFind.findRepresentative(posnA), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnB), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnC), posnC);
    t.checkExpect(posnUnionFind.findRepresentative(posnD), posnD);

    // Connect posnC to posnD and then that tree to posnA
    posnUnionFind.union(posnC, posnD);
    posnUnionFind.union(posnA, posnC);

    // Check the new representatives
    t.checkExpect(posnUnionFind.findRepresentative(posnA), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnB), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnC), posnA);
    t.checkExpect(posnUnionFind.findRepresentative(posnD), posnA);
  }

  void testAddFindAddBehavior(Tester t) {
    this.initTestData();

    AddFind<String> stringAddFind = new AddFind<String>();

    // Add a string as a root
    stringAddFind.addRoot("ROOT");

    // It is now contained in the AddFind
    t.checkExpect(stringAddFind.contains("ROOT"), true);

    // Show that attempting to add another value as the root fails
    t.checkException(new UnsupportedOperationException("Cannot add a root to a non-empty AddFind"),
        stringAddFind, "addRoot", "ROOT II");

    // Add a value to the AddFind and show that
    // the root still remains the root
    stringAddFind.add("STEM1", "ROOT");

    // It is now contained in the AddFind (as well as "ROOT")
    t.checkExpect(stringAddFind.contains("ROOT"), true);
    t.checkExpect(stringAddFind.contains("STEM1"), true);

    // Add some more values and check the tree we hare building
    // up. We are building the following in-tree:
    //
    // ROOT is at the ROOT
    // STEM1 branches from ROOT
    // STEM2 branches from ROOT
    // STEM3 branches from ROOT
    //
    // S1STEM1 branches from STEM1
    // S1STEM2 branches from STEM1
    // S2STEM1 branches from STEM2
    // S2STEM2 branches from STEM2
    // S2S2STEM1 branches form S2STEM2
    stringAddFind.add("STEM2", "ROOT");
    stringAddFind.add("STEM3", "ROOT");

    stringAddFind.add("S1STEM1", "STEM1");
    stringAddFind.add("S1STEM2", "STEM1");
    stringAddFind.add("S2STEM1", "STEM2");
    stringAddFind.add("S2STEM2", "STEM2");

    stringAddFind.add("S2S2STEM1", "S2STEM2");

    // Trying to add an element to itself in the
    // AddFind is not allowed: elements must connect to nodes
    // already in the in-tree
    t.checkException(
        new IllegalArgumentException("Cannot directly add an item to the `AddFind<T>` "
            + "whose parent is itself since that element is not the root"),
        stringAddFind, "add", "S2S2STEM1", "S2S2STEM1");

    t.checkException(
        new IllegalArgumentException("Cannot directly add an item to the `AddFind<T>` "
            + "whose parent is itself since that element is not the root"),
        stringAddFind, "add", "STEMSAME", "STEMSAME");

    // Trying to add an element that is already a part of the AddFind is an
    // exception
    t.checkException(
        new IllegalArgumentException(
            "Cannot add the given element to the AddFind because it is already in the AddFind"),
        stringAddFind, "add", "S2STEM1", "STEM1");
    t.checkException(
        new IllegalArgumentException(
            "Cannot add the given element to the AddFind because it is already in the AddFind"),
        stringAddFind, "add", "S2STEM2", "STEM1");

    // Trying to add something to a parent that is not already in the AddFind
    // is an exception
    t.checkException(new IllegalArgumentException(
        "Cannot add an element to the AddFind to a node that is not already part of the AddFind"),
        stringAddFind, "add", "S2S2STEM2", "S2STEM3");

    t.checkException(new IllegalArgumentException(
        "Cannot add an element to the AddFind to a node that is not already part of the AddFind"),
        stringAddFind, "add", "S2S2STEM2", "S2STEM4");

    // Check everything is in the AddFind
    t.checkExpect(stringAddFind.contains("ROOT"), true);
    t.checkExpect(stringAddFind.contains("STEM1"), true);
    t.checkExpect(stringAddFind.contains("STEM2"), true);
    t.checkExpect(stringAddFind.contains("STEM3"), true);
    t.checkExpect(stringAddFind.contains("S1STEM1"), true);
    t.checkExpect(stringAddFind.contains("S1STEM2"), true);
    t.checkExpect(stringAddFind.contains("S2STEM1"), true);
    t.checkExpect(stringAddFind.contains("S2STEM2"), true);
    t.checkExpect(stringAddFind.contains("S2S2STEM1"), true);
  }

  void testAddFindGetLineage(Tester t) {
    AddFind<String> stringAddFind = new AddFind<String>();

    // We are building the following in-tree (same as the previous method)
    //
    // ROOT is at the ROOT
    // STEM1 branches from ROOT
    // STEM2 branches from ROOT
    // STEM3 branches from ROOT
    //
    // S1STEM1 branches from STEM1
    // S1STEM2 branches from STEM1
    // S2STEM1 branches from STEM2
    // S2STEM2 branches from STEM2
    // S2S2STEM1 branches form S2STEM2
    stringAddFind.addRoot("ROOT");
    stringAddFind.add("STEM1", "ROOT");
    stringAddFind.add("STEM2", "ROOT");
    stringAddFind.add("STEM3", "ROOT");

    stringAddFind.add("S1STEM1", "STEM1");
    stringAddFind.add("S1STEM2", "STEM1");
    stringAddFind.add("S2STEM1", "STEM2");
    stringAddFind.add("S2STEM2", "STEM2");
    stringAddFind.add("S2S2STEM1", "S2STEM2");

    // Now, test the pathways up the AddFind
    t.checkExpect(stringAddFind.getLineage("ROOT"), new ArrayList<String>(List.of("ROOT")));
    t.checkExpect(stringAddFind.getLineage("STEM1"),
        new ArrayList<String>(List.of("STEM1", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("STEM2"),
        new ArrayList<String>(List.of("STEM2", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("STEM3"),
        new ArrayList<String>(List.of("STEM3", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("S1STEM1"),
        new ArrayList<String>(List.of("S1STEM1", "STEM1", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("S1STEM2"),
        new ArrayList<String>(List.of("S1STEM2", "STEM1", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("S2STEM1"),
        new ArrayList<String>(List.of("S2STEM1", "STEM2", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("S2STEM2"),
        new ArrayList<String>(List.of("S2STEM2", "STEM2", "ROOT")));
    t.checkExpect(stringAddFind.getLineage("S2S2STEM1"),
        new ArrayList<String>(List.of("S2S2STEM1", "S2STEM2", "STEM2", "ROOT")));

    // Attempting to get the lineage of a node that is not
    // in the AddFind is an exception
    t.checkException(
        new IllegalArgumentException("The given element is not present in the AddFind"),
        stringAddFind, "getLineage", "NOT_IN_TREE");
    t.checkException(
        new IllegalArgumentException("The given element is not present in the AddFind"),
        stringAddFind, "getLineage", "NOT_IN_TREE_AGAIN");

    // Try another AddFind of a different type: Posns
    AddFind<Posn> posnAddFind = new AddFind<Posn>();

    Posn posA = new Posn(10, 10);
    Posn posAA = new Posn(11, 1);
    Posn posAB = new Posn(5, 1);
    Posn posAC = new Posn(6, -9);
    Posn posACA = new Posn(3, 5);
    Posn posACB = new Posn(19, -20);

    // Connect them together
    posnAddFind.addRoot(posA);
    posnAddFind.add(posAA, posA);
    posnAddFind.add(posAB, posA);
    posnAddFind.add(posAC, posA);
    posnAddFind.add(posACA, posAC);
    posnAddFind.add(posACB, posAC);

    // Test all of the paths through the AddFind
    t.checkExpect(posnAddFind.getLineage(posA), new ArrayList<Posn>(List.of(posA)));
    t.checkExpect(posnAddFind.getLineage(posAA), new ArrayList<Posn>(List.of(posAA, posA)));
    t.checkExpect(posnAddFind.getLineage(posAB), new ArrayList<Posn>(List.of(posAB, posA)));
    t.checkExpect(posnAddFind.getLineage(posAC), new ArrayList<Posn>(List.of(posAC, posA)));
    t.checkExpect(posnAddFind.getLineage(posACA),
        new ArrayList<Posn>(List.of(posACA, posAC, posA)));
    t.checkExpect(posnAddFind.getLineage(posACB),
        new ArrayList<Posn>(List.of(posACB, posAC, posA)));
  }

  void testAddFindDepthOf(Tester t) {
    this.initTestData();

    AddFind<String> stringAddFind = new AddFind<String>();

    // We are building the following in-tree (same as before)
    //
    // ROOT is at the ROOT
    // STEM1 branches from ROOT
    // STEM2 branches from ROOT
    // STEM3 branches from ROOT
    //
    // S1STEM1 branches from STEM1
    // S1STEM2 branches from STEM1
    // S2STEM1 branches from STEM2
    // S2STEM2 branches from STEM2
    // S2S2STEM1 branches form S2STEM2
    stringAddFind.addRoot("ROOT");
    stringAddFind.add("STEM1", "ROOT");
    stringAddFind.add("STEM2", "ROOT");
    stringAddFind.add("STEM3", "ROOT");

    stringAddFind.add("S1STEM1", "STEM1");
    stringAddFind.add("S1STEM2", "STEM1");
    stringAddFind.add("S2STEM1", "STEM2");
    stringAddFind.add("S2STEM2", "STEM2");

    stringAddFind.add("S2S2STEM1", "S2STEM2");

    t.checkExpect(stringAddFind.depthOf("ROOT"), 0);
    t.checkExpect(stringAddFind.depthOf("STEM1"), 1);
    t.checkExpect(stringAddFind.depthOf("STEM2"), 1);
    t.checkExpect(stringAddFind.depthOf("STEM3"), 1);
    t.checkExpect(stringAddFind.depthOf("S1STEM1"), 2);
    t.checkExpect(stringAddFind.depthOf("S1STEM2"), 2);
    t.checkExpect(stringAddFind.depthOf("S2STEM1"), 2);
    t.checkExpect(stringAddFind.depthOf("S2STEM2"), 2);
    t.checkExpect(stringAddFind.depthOf("S2S2STEM1"), 3);

    // Attempting to compute the depth of a node not in the tree is an exception
    t.checkException(
        new IllegalArgumentException("The given element is not in the AddFind to get a depth for"),
        stringAddFind, "depthOf", "NOT_IN_TREE");
    t.checkException(
        new IllegalArgumentException("The given element is not in the AddFind to get a depth for"),
        stringAddFind, "depthOf", "NOT_IN_TREE_AGAIN");
  }

  void testAddFindIsEmpty(Tester t) {
    this.initTestData();

    AddFind<String> stringAddFind = new AddFind<String>();

    // Nothing is in the AddFind
    t.checkExpect(stringAddFind.isEmpty(), true);

    // Add something as the root, and now we have something
    stringAddFind.addRoot("A");
    t.checkExpect(stringAddFind.isEmpty(), false);
  }

  void testAddFindIterator(Tester t) {
    this.initTestData();

    // The order in which we receive the
    // items in an AddFind as we iterate over its structure
    // is not guaranteed (as documented with its purpose
    // statement). However, we can at least check if we've
    // hit everything in the AddFind
    AddFind<String> addFindString = new AddFind<String>();

    // If we get an iterator right now, it won't be very useful
    Iterator<String> notUseful = addFindString.iterator();

    t.checkExpect(notUseful.hasNext(), false);

    // If we add a single item, we expect to go over only that item
    addFindString.addRoot("ONE_ITEM");
    Iterator<String> oneItemIter = addFindString.iterator();
    t.checkExpect(oneItemIter.hasNext(), true);
    t.checkExpect(oneItemIter.next(), "ONE_ITEM");
    t.checkExpect(oneItemIter.hasNext(), false);

    // If we add more items, things get a little trickier. The order
    // is not guaranteed. Here, we just check that we've seen both
    // of the items
    addFindString.add("ANOTHER_ITEM", "ONE_ITEM");
    Iterator<String> twoItemIter = addFindString.iterator();

    ArrayList<String> seenItems = new ArrayList<String>();
    while (twoItemIter.hasNext()) {
      seenItems.add(twoItemIter.next());
    }

    t.checkExpect(seenItems.contains("ONE_ITEM"), true);
    t.checkExpect(seenItems.contains("ANOTHER_ITEM"), true);

    // Add some more
    addFindString.add("THIRD_ITEM", "ONE_ITEM");
    addFindString.add("FOURTH_ITEM", "THIRD_ITEM");
    addFindString.add("FIFTH_ITEM", "THIRD_ITEM");
    addFindString.add("SIXTH_ITEM", "FIFTH_ITEM");
    Iterator<String> manyItemIter = addFindString.iterator();

    seenItems = new ArrayList<String>();
    while (manyItemIter.hasNext()) {
      seenItems.add(manyItemIter.next());
    }
    t.checkExpect(seenItems.contains("ONE_ITEM"), true);
    t.checkExpect(seenItems.contains("ANOTHER_ITEM"), true);
    t.checkExpect(seenItems.contains("THIRD_ITEM"), true);
    t.checkExpect(seenItems.contains("FOURTH_ITEM"), true);
    t.checkExpect(seenItems.contains("FIFTH_ITEM"), true);
    t.checkExpect(seenItems.contains("SIXTH_ITEM"), true);
  }

  void testSearchOneStep(Tester t) {

    // using a self made grid that looks like this

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    // we'll make three snapshots, each for each type of search

    ASearchSnapshot snapManual = new ManualSearchSnapshot(cell1, cell9);
    ASearchSnapshot snapAuto = new DepthFirstSnapshot(cell1, cell9);

    // we'll start by searching manually. Manual searchOneStep only returns a
    // boolean, with no side effects, so we'll see if it returns true when it
    // reaches the finish

    snapManual.maybeMoveDown();
    t.checkExpect(snapManual.searchOneStep(), false);
    snapManual.maybeMoveDown();
    t.checkExpect(snapManual.searchOneStep(), false);
    snapManual.maybeMoveRight();
    t.checkExpect(snapManual.searchOneStep(), false);
    snapManual.maybeMoveRight();
    t.checkExpect(snapManual.searchOneStep(), true);

    // for manual searches, our current place will stay the same even after we
    // finish, assuming we don't move further. So we can continue calling
    // search one step, returning true each time.
    t.checkExpect(snapManual.searchOneStep(), true);

    // for depth first, we should cannot manually move our current cell,
    // as searchOneStep has the side effect of searching one step for us.
    // We should find the finish after 5 steps.

    t.checkExpect(snapAuto.searchOneStep(), false);
    t.checkExpect(snapAuto.searchOneStep(), false);
    t.checkExpect(snapAuto.searchOneStep(), false);
    t.checkExpect(snapAuto.searchOneStep(), false);
    t.checkExpect(snapAuto.searchOneStep(), true);

    // after we finish, we haven't seen everything, so we can still search.
    // SearchOneStep will return false, however, as we are no longer at our finish,
    // so this kind of usage should not occur
    t.checkExpect(snapAuto.searchOneStep(), false);

  }

  void testHasMoreToSee(Tester t) {

    // same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    ASearchSnapshot snapManual = new ManualSearchSnapshot(cell1, cell9);
    ASearchSnapshot snapAuto = new DepthFirstSnapshot(cell1, cell9);

    // for manual snapshots, we should always be able to see more, so now amount
    // of travel can eliminate what we have to see
    t.checkExpect(snapManual.hasMoreToSee(), true);
    snapManual.maybeMoveDown();
    snapManual.maybeMoveDown();
    t.checkExpect(snapManual.hasMoreToSee(), true);
    snapManual.maybeMoveRight();
    snapManual.maybeMoveRight();
    t.checkExpect(snapManual.hasMoreToSee(), true);
    snapManual.maybeMoveLeft();
    snapManual.maybeMoveLeft();
    t.checkExpect(snapManual.hasMoreToSee(), true);
    snapManual.maybeMoveUp();
    snapManual.maybeMoveUp();
    t.checkExpect(snapManual.hasMoreToSee(), true);

    // for automatic snapshots, we also should not run out of things to see,
    // as we are guaranteed to find the end before we run out of items,
    // as long as we correctly processed our minimum spanning tree

    t.checkExpect(snapAuto.hasMoreToSee(), true);
    // we know from previous tests that depth first solves this maze
    // in 5 steps
    snapAuto.searchOneStep();
    snapAuto.searchOneStep();
    snapAuto.searchOneStep();
    snapAuto.searchOneStep();
    snapAuto.searchOneStep();
    // we did not run out by the time we finished
    t.checkExpect(snapAuto.hasMoreToSee(), true);

    // while we cannot actually search further, we can see that attempting
    // a continued search does not leave us unable to see more, even after
    // we have finished searching
    snapAuto.searchOneStep();
    t.checkExpect(snapAuto.hasMoreToSee(), true);

  }

  void testDistanceToFathest(Tester t) {

    // using the same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    ASearchSnapshot snapDepth = new DepthFirstSnapshot(cell1, cell9);
    ASearchSnapshot snapBreadth = new BreadthFirstSnapshot(cell1, cell9);

    // after a two steps, we should have a farthest distance of one
    snapDepth.searchOneStep();
    snapDepth.searchOneStep();
    t.checkExpect(snapDepth.distanceToFarthest(), 1);
    // and another
    snapDepth.searchOneStep();
    t.checkExpect(snapDepth.distanceToFarthest(), 2);
    // and another
    snapDepth.searchOneStep();
    t.checkExpect(snapDepth.distanceToFarthest(), 3);
    // and another
    snapDepth.searchOneStep();
    t.checkExpect(snapDepth.distanceToFarthest(), 4);
    // the only cell in the grid that is 4 away from the start is the finish
    // cell, which works out well, as we expect the depth first search to
    // travel directly to the target cell, given our maze design

    // breadth first operates a bit differently:
    // we begin as normal
    snapBreadth.searchOneStep();
    snapBreadth.searchOneStep();
    t.checkExpect(snapBreadth.distanceToFarthest(), 1);
    // but at this step we travel in a different direction, so the
    // farthest distance stays the same
    snapBreadth.searchOneStep();
    t.checkExpect(snapBreadth.distanceToFarthest(), 1);
    // now that we've exhausted the start connections, we'll move to
    // cells that are 2 away
    snapBreadth.searchOneStep();
    t.checkExpect(snapBreadth.distanceToFarthest(), 2);
    snapBreadth.searchOneStep();
    t.checkExpect(snapBreadth.distanceToFarthest(), 2);
    snapBreadth.searchOneStep();
    t.checkExpect(snapBreadth.distanceToFarthest(), 2);
    // this continues for the rest of the maze

  }

  void testMakeAllSeenAndUnseen(Tester t) {

    // covers makeAllSeen and makeAllUnseen for an automatic search
    // using the same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    // this snapshot takes 5 moves to reach the end
    ASearchSnapshot snap = new DepthFirstSnapshot(cell1, cell9);
    snap.searchOneStep();
    snap.searchOneStep();
    snap.searchOneStep();
    snap.searchOneStep();
    t.checkExpect(snap.searchOneStep(), true);

    // now we make all the one's we've passed through seen
    snap.makeAllSeen();

    t.checkExpect(cell1.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, true, false, false));
    t.checkExpect(cell2.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, false));
    t.checkExpect(cell3.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, true, false));
    t.checkExpect(cell4.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, false));
    t.checkExpect(cell5.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, true));
    t.checkExpect(cell7.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, false, true, false, true));
    t.checkExpect(cell8.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, true));

    snap.makeAllUnseen();

    t.checkExpect(cell1.toImage(10), this.boxedImage2(10, Constants.CELL_START_COLOR,
        Constants.WALL_COLOR, true, true, false, false));
    t.checkExpect(cell2.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, false));
    t.checkExpect(cell3.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, true, false));
    t.checkExpect(cell4.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, false));
    t.checkExpect(cell5.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, true));
    t.checkExpect(cell7.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, false, true));
    t.checkExpect(cell8.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, true));
  }

  void testMakeAllSeenAndUnseenPt2(Tester t) {

    // same as the last method, but on a manual search, so we can see some
    // items that aren't in our path

    // using the same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    // this snapshot takes 5 moves to reach the end
    ASearchSnapshot snap = new ManualSearchSnapshot(cell1, cell9);
    snap.maybeMoveRight();
    snap.maybeMoveRight();
    snap.maybeMoveLeft();
    snap.maybeMoveLeft();
    snap.maybeMoveDown();
    snap.maybeMoveDown();
    snap.maybeMoveRight();
    snap.maybeMoveRight();

    // now we make all the one's we've passed through seen
    snap.makeAllSeen();

    t.checkExpect(cell1.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, true, false, false));
    t.checkExpect(cell2.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, false));
    t.checkExpect(cell3.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, true, false));
    t.checkExpect(cell4.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, false));
    t.checkExpect(cell5.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, true));
    t.checkExpect(cell7.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, false, true, false, true));
    t.checkExpect(cell8.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, true));

    snap.makeAllUnseen();

    t.checkExpect(cell1.toImage(10), this.boxedImage2(10, Constants.CELL_START_COLOR,
        Constants.WALL_COLOR, true, true, false, false));
    t.checkExpect(cell2.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, false));
    t.checkExpect(cell3.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, true, false));
    t.checkExpect(cell4.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, false));
    t.checkExpect(cell5.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, true));
    t.checkExpect(cell7.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, false, true));
    t.checkExpect(cell8.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, true));
  }

  void testMakeAllInPath(Tester t) {

    // using the same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    ASearchSnapshot snapManual = new ManualSearchSnapshot(cell1, cell9);

    // we'll complete a manual search with some bad moves, to be able to
    // differentiate between unseen, seen, and in path in a single method.
    snapManual.maybeMoveRight();
    snapManual.maybeMoveRight();
    snapManual.maybeMoveLeft();
    snapManual.maybeMoveLeft();
    snapManual.maybeMoveDown();
    snapManual.maybeMoveDown();
    snapManual.maybeMoveRight();
    snapManual.maybeMoveRight();
    t.checkExpect(snapManual.searchOneStep(), true);
    snapManual.makeAllSeen();
    snapManual.makeAllInPath();

    t.checkExpect(cell1.toImage(10), this.boxedImage2(10, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR, true, true, false, false));
    t.checkExpect(cell2.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, false, false));
    t.checkExpect(cell3.toImage(10), this.boxedImage2(10, Constants.CELL_SEEN_COLOR,
        Constants.WALL_COLOR, true, false, true, false));
    t.checkExpect(cell4.toImage(10), this.boxedImage2(10, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR, false, true, true, false));
    t.checkExpect(cell5.toImage(10), this.boxedImage2(10, Constants.CELL_UNSEEN_COLOR,
        Constants.WALL_COLOR, false, true, true, true));
    t.checkExpect(cell7.toImage(10), this.boxedImage2(10, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR, false, true, false, true));
    t.checkExpect(cell8.toImage(10), this.boxedImage2(10, Constants.CELL_PATH_COLOR,
        Constants.WALL_COLOR, true, false, false, true));

  }

  void testDistanceFromStart(Tester t) {

    // using the same grid as before

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    ASearchSnapshot snap = new DepthFirstSnapshot(cell1, cell9);
    snap.searchOneStep();
    snap.searchOneStep();
    snap.searchOneStep();
    snap.searchOneStep();
    t.checkExpect(snap.searchOneStep(), true);

    // check all of the items in our true path
    t.checkExpect(snap.distanceFromStart(cell1), 0);
    t.checkExpect(snap.distanceFromStart(cell4), 1);
    t.checkExpect(snap.distanceFromStart(cell7), 2);
    t.checkExpect(snap.distanceFromStart(cell8), 3);
    t.checkExpect(snap.distanceFromStart(cell9), 4);

    // check for exceptions on some cells that we haven't seen
    t.checkException(
        new IllegalArgumentException(
            "The search algorithm state has not yet processed the given ICell"),
        snap, "distanceFromStart", cell2);
    t.checkException(
        new IllegalArgumentException(
            "The search algorithm state has not yet processed the given ICell"),
        snap, "distanceFromStart", cell5);
    t.checkException(
        new IllegalArgumentException(
            "The search algorithm state has not yet processed the given ICell"),
        snap, "distanceFromStart", cell6);

  }

  void testMaybeMoveAllDirections(Tester t) {

    // covers tests for maybeMoveUp, Down, Left, and Right

    // we will be using this self made grid, which looks like this

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    // while we cannot actually see where we are, we are able to tell
    // when we have reached our finish line. We will try a few different
    // paths, until we reach the finish line

    // search 1: direct path

    ASearchSnapshot snap = new ManualSearchSnapshot(cell1, cell9);
    t.checkExpect(snap.searchOneStep(), false);
    snap.maybeMoveDown();
    t.checkExpect(snap.searchOneStep(), false);
    snap.maybeMoveDown();
    t.checkExpect(snap.searchOneStep(), false);
    snap.maybeMoveRight();
    t.checkExpect(snap.searchOneStep(), false);
    snap.maybeMoveRight();
    t.checkExpect(snap.searchOneStep(), true);

    // search 2: direct path, with a number of illegal moves
    ASearchSnapshot snap2 = new ManualSearchSnapshot(cell1, cell9);
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveUp();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveLeft();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveDown();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveRight();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveDown();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveRight();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveUp();
    t.checkExpect(snap2.searchOneStep(), false);
    snap2.maybeMoveRight();
    t.checkExpect(snap2.searchOneStep(), true);
  }

  void testAllOrAdjacentOutCells(Tester t) {

    // tests for AllAdjacentCells and AdjacentOutCells
    // we will use the same grid from the previous test

    /*   _____
     *  |1 2 3|
     *  |3|_|_| <-- 5, 6
     *  |7 8 9|
     *   -----
     *
     */

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {
      edge.becomeEdgeInTree();
    }

    t.checkExpect(cell5.allAdjacentCells(),
        new ArrayList<ICell>(List.of(cell2, cell4, cell6, cell8)));

    // this grid is not big enough to test further all adjacent edges,
    // but more tests for allAdjacentCells() can be found in each method
    // for testing connections, where the method is tested alongside the
    // connection methods

    t.checkExpect(cell1.adjacentOutCells(), new ArrayList<ICell>(List.of(cell2, cell4)));
    t.checkExpect(cell2.adjacentOutCells(), new ArrayList<ICell>(List.of(cell1, cell3, cell5)));
    t.checkExpect(cell3.adjacentOutCells(), new ArrayList<ICell>(List.of(cell2, cell6)));
    t.checkExpect(cell4.adjacentOutCells(), new ArrayList<ICell>(List.of(cell1, cell7)));
    t.checkExpect(cell5.adjacentOutCells(), new ArrayList<ICell>(List.of(cell2)));
    t.checkExpect(cell7.adjacentOutCells(), new ArrayList<ICell>(List.of(cell4, cell8)));
    t.checkExpect(cell8.adjacentOutCells(), new ArrayList<ICell>(List.of(cell7, cell9)));
    t.checkExpect(cell9.adjacentOutCells(), new ArrayList<ICell>(List.of(cell8)));
  }

  void testMaybeAllDirections(Tester t) {

    // includes tests for maybeTop, maybeBottom, maybeLeft, and maybeRight

    Cell cell1 = new Cell();
    Cell cell2 = new Cell();
    Cell cell3 = new Cell();
    Cell cell4 = new Cell();
    Cell cell5 = new Cell();
    Cell cell6 = new Cell();
    Cell cell7 = new Cell();
    Cell cell8 = new Cell();
    Cell cell9 = new Cell();

    // connect all cells horizontally
    CellEdge edge12 = cell1.connectAtRight(cell2, 4);
    CellEdge edge23 = cell2.connectAtRight(cell3, 4);
    CellEdge edge45 = cell4.connectAtRight(cell5, 100);
    CellEdge edge56 = cell5.connectAtRight(cell6, 100);
    CellEdge edge78 = cell7.connectAtRight(cell8, 4);
    CellEdge edge89 = cell8.connectAtRight(cell9, 4);

    // connect all cells vertically
    CellEdge edge14 = cell1.connectAtBottom(cell4, 4);
    CellEdge edge25 = cell2.connectAtBottom(cell5, 4);
    CellEdge edge36 = cell3.connectAtBottom(cell6, 4);
    CellEdge edge47 = cell4.connectAtBottom(cell7, 4);
    CellEdge edge58 = cell5.connectAtBottom(cell8, 100);
    CellEdge edge69 = cell6.connectAtBottom(cell6, 100);

    // we'll then assemble lists for allPathways
    ArrayList<CellEdge> edges = new ArrayList<CellEdge>(List.of(edge12, edge23, edge45, edge56,
        edge78, edge89, edge14, edge25, edge36, edge47, edge58, edge69));

    ArrayList<ICell> cells = new ArrayList<ICell>(
        List.of(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9));

    ArrayList<CellEdge> openEdges = new MazeUtils().allPathways(edges, cells);

    // go over all edges that are open and set them to open
    // Note: we cannot use breakAllWalls here, because we do not
    // technically have a maze, only a maze-looking grid
    for (CellEdge edge : openEdges) {

      edge.becomeEdgeInTree();
    }

    // for blocked paths, we return the cell that requested its neighbor
    t.checkExpect(cell1.maybeTop(), cell1);
    t.checkExpect(cell1.maybeLeft(), cell1);
    t.checkExpect(cell5.maybeBottom(), cell5);
    t.checkExpect(cell5.maybeRight(), cell5);
    t.checkExpect(cell6.maybeRight(), cell6);
    t.checkExpect(cell8.maybeTop(), cell8);
    t.checkExpect(cell9.maybeBottom(), cell9);

    // for open paths, we return the appropriate neighbor
    t.checkExpect(cell1.maybeBottom(), cell4);
    t.checkExpect(cell1.maybeRight(), cell2);
    t.checkExpect(cell8.maybeRight(), cell9);
    t.checkExpect(cell9.maybeLeft(), cell8);
    t.checkExpect(cell4.maybeBottom(), cell7);
    t.checkExpect(cell4.maybeTop(), cell1);
    t.checkExpect(cell3.maybeLeft(), cell2);
  }

  void testIsOpenPath(Tester t) {
    t.checkExpect(new OpenState().isOpenPath(), true);
    t.checkExpect(new ClosedState().isOpenPath(), false);
  }

  void testNextHorizontalWeight(Tester t) {
    this.initTestData();
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, 11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", 10, -11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, -11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", 10, -11);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, 11);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, -11);
    EdgeWeightGen weightGen = new EdgeWeightGen(100, 100, 0, 0);

    t.checkExpect(weightGen.nextHorizontalWeight(), 67);
    t.checkExpect(weightGen.nextHorizontalWeight(), 72);
    t.checkExpect(weightGen.nextHorizontalWeight(), 93);
    t.checkExpect(weightGen.nextHorizontalWeight(), 5);
    t.checkExpect(weightGen.nextHorizontalWeight(), 9);
    t.checkExpect(weightGen.nextHorizontalWeight(), 54);
    t.checkExpect(weightGen.nextHorizontalWeight(), 82);
    t.checkExpect(weightGen.nextHorizontalWeight(), 42);
    // Try a different seed and different maximum value
    weightGen = new EdgeWeightGen(400, 200, 20, 0);
    t.checkExpect(weightGen.nextHorizontalWeight(), 330);
    t.checkExpect(weightGen.nextHorizontalWeight(), 287);
    t.checkExpect(weightGen.nextHorizontalWeight(), 286);
    t.checkExpect(weightGen.nextHorizontalWeight(), 202);
    t.checkExpect(weightGen.nextHorizontalWeight(), 293);
    t.checkExpect(weightGen.nextHorizontalWeight(), 239);
    t.checkExpect(weightGen.nextHorizontalWeight(), 106);
    t.checkExpect(weightGen.nextHorizontalWeight(), 152);
  }

  void testNextVerticalWeight(Tester t) {
    this.initTestData();
    this.initTestData();
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, 11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", 10, -11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, -11, 19120, 90);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", 10, -11);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, 11);
    t.checkConstructorException(
        new IllegalArgumentException("Maximum weights must be greater than or equal to zero."),
        "EdgeWeightGen", -10, -11);
    EdgeWeightGen weightGen = new EdgeWeightGen(100, 100, 0, 0);
    t.checkExpect(weightGen.nextVerticalWeight(), 67);
    t.checkExpect(weightGen.nextVerticalWeight(), 72);
    t.checkExpect(weightGen.nextVerticalWeight(), 93);
    t.checkExpect(weightGen.nextVerticalWeight(), 5);
    t.checkExpect(weightGen.nextVerticalWeight(), 9);
    t.checkExpect(weightGen.nextVerticalWeight(), 54);
    t.checkExpect(weightGen.nextVerticalWeight(), 82);
    t.checkExpect(weightGen.nextVerticalWeight(), 42);
    // Try a different seed and different maximum value
    weightGen = new EdgeWeightGen(400, 200, 20, 0);
    t.checkExpect(weightGen.nextVerticalWeight(), 102);
    t.checkExpect(weightGen.nextVerticalWeight(), 34);
    t.checkExpect(weightGen.nextVerticalWeight(), 139);
    t.checkExpect(weightGen.nextVerticalWeight(), 56);
    t.checkExpect(weightGen.nextVerticalWeight(), 149);
    t.checkExpect(weightGen.nextVerticalWeight(), 158);
    t.checkExpect(weightGen.nextVerticalWeight(), 23);
    t.checkExpect(weightGen.nextVerticalWeight(), 63);
  }

  void testClamp(Tester t) {
    this.initTestData();

    t.checkInexact(new DoubleUtils().clamp(10.0, 11.0, 12.0), 11.0, 0.001);
    t.checkInexact(new DoubleUtils().clamp(-100.0, -6.0, 5.0), -6.0, 0.001);
    t.checkInexact(new DoubleUtils().clamp(8.7, 4.0, 12.0), 8.7, 0.001);
    t.checkInexact(new DoubleUtils().clamp(14.0, 0.0, 12.0), 12.0, 0.001);
    t.checkInexact(new DoubleUtils().clamp(-10.0, -12.0, -10.0), -10.0, 0.001);

    t.checkException(
        new IllegalArgumentException(
            "Cannot clamp a value in a range where the lower value is greater than the higher one"),
        new DoubleUtils(), "clamp", -10.0, -13.0, -15.0);
    t.checkException(
        new IllegalArgumentException(
            "Cannot clamp a value in a range where the lower value is greater than the higher one"),
        new DoubleUtils(), "clamp", -10.0, 3.0, 2.0);
    t.checkException(
        new IllegalArgumentException(
            "Cannot clamp a value in a range where the lower value is greater than the higher one"),
        new DoubleUtils(), "clamp", -9.0, 11.0, 2.0);
  }

  void testColorUtilsInterpolate(Tester t) {
    this.initTestData();

    t.checkExpect(new ColorUtils().interpolate(-0.4, Color.RED, Color.BLUE), new Color(255, 0, 0));
    t.checkExpect(new ColorUtils().interpolate(0.0, Color.RED, Color.BLUE), new Color(255, 0, 0));
    t.checkExpect(new ColorUtils().interpolate(0.2, Color.RED, Color.BLUE), new Color(204, 0, 51));
    t.checkExpect(new ColorUtils().interpolate(0.5, Color.RED, Color.BLUE), new Color(127, 0, 127));
    t.checkExpect(new ColorUtils().interpolate(0.8, Color.RED, Color.BLUE), new Color(50, 0, 204));
    t.checkExpect(new ColorUtils().interpolate(0.9, Color.RED, Color.BLUE), new Color(25, 0, 229));
    t.checkExpect(new ColorUtils().interpolate(1.0, Color.RED, Color.BLUE), new Color(0, 0, 255));
    t.checkExpect(new ColorUtils().interpolate(1.2, Color.RED, Color.BLUE), new Color(0, 0, 255));
  }

  void testMap(Tester t) {
    this.initTestData();

    t.checkExpect(new ArrayUtils().map(new ArrayList<String>(), new ToLength()),
        new ArrayList<Integer>());
    t.checkExpect(new ArrayUtils().map(new ArrayList<String>(List.of("A", "B")), new ToLength()),
        new ArrayList<Integer>(List.of(1, 1)));
    t.checkExpect(
        new ArrayUtils().map(new ArrayList<String>(List.of("9876", "12345", "^")), new ToLength()),
        new ArrayList<Integer>(List.of(4, 5, 1)));
    t.checkExpect(new ArrayUtils().map(new ArrayList<Integer>(List.of(-10, 11, 10)),
        new IntegerAdditiveInverse()), new ArrayList<Integer>(List.of(10, -11, -10)));
  }

  void testFilter(Tester t) {
    this.initTestData();

    t.checkExpect(
        new ArrayUtils().filter(new ArrayList<Integer>(List.of(-10, 11, 10)), new IsPositive()),
        new ArrayList<Integer>(List.of(11, 10)));

    t.checkExpect(new ArrayUtils().filter(new ArrayList<Integer>(List.of(-10, -11, -12, -10, 0)),
        new IsPositive()), new ArrayList<Integer>());

    t.checkExpect(
        new ArrayUtils().filter(new ArrayList<Integer>(List.of(-10, 11, 10)), new IsPositive()),
        new ArrayList<Integer>(List.of(11, 10)));

    t.checkExpect(new ArrayUtils().filter(new ArrayList<Integer>(), new IsPositive()),
        new ArrayList<Integer>());

    t.checkExpect(
        new ArrayUtils().filter(
            new ArrayList<Posn>(List.of(new Posn(-100, 10), new Posn(100, 10), new Posn(10, 10))),
            new InFirstQuadrant()),
        new ArrayList<Posn>(List.of(new Posn(100, 10), new Posn(10, 10))));
  }

  void testFoldl(Tester t) {
    this.initTestData();

    t.checkExpect(
        new ArrayUtils().foldl(new ArrayList<String>(List.of("A", "B")), "D", new CombineStrings()),
        "DAB");
    t.checkExpect(new ArrayUtils().foldl(new ArrayList<String>(List.of("E", "M", "E")), "M",
        new CombineStrings()), "MEME");
    t.checkExpect(new ArrayUtils().foldl(new ArrayList<String>(), "base", new CombineStrings()),
        "base");
  }

  void testFirstToSatisfy(Tester t) {
    this.initTestData();

    t.checkExpect(new ArrayUtils().firstToSatisfy(new ArrayList<Integer>(List.of(-10, 11, 10)),
        new IsPositive()), 11);

    t.checkExpect(new ArrayUtils()
        .firstToSatisfy(new ArrayList<Integer>(List.of(-10, -10, -100, 12)), new IsPositive()), 12);

    t.checkExpect(new ArrayUtils().firstToSatisfy(new ArrayList<Integer>(List.of(-10, 11, 10)),
        new IsNegativeOrZero()), -10);

    t.checkException(new IllegalArgumentException("Cannot find an item in an empty list"),
        new ArrayUtils(), "firstToSatisfy", new ArrayList<Integer>(), new IsPositive());

    t.checkException(
        new IllegalArgumentException("No argument in the list satisfied the given" + " function."),
        new ArrayUtils(), "firstToSatisfy", new ArrayList<Integer>(List.of(-10, -11, -100)),
        new IsPositive());
  }

  void testICellStateStandardColor(Tester t) {
    this.initTestData();
    t.checkExpect(new UnseenState().standardColor(), Constants.CELL_UNSEEN_COLOR);
    t.checkExpect(new SeenState().standardColor(), Constants.CELL_SEEN_COLOR);
    t.checkExpect(new StartState().standardColor(), Constants.CELL_START_COLOR);
    t.checkExpect(new FinishState().standardColor(), Constants.CELL_FINISH_COLOR);
    t.checkExpect(new PathState().standardColor(), Constants.CELL_PATH_COLOR);
  }

  void testICellStateInterpolatedColor(Tester t) {
    this.initTestData();

    t.checkExpect(new UnseenState().interpolatedColor(10, 100),
        new ColorUtils().interpolate(0.1, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new UnseenState().interpolatedColor(-10, 100),
        new ColorUtils().interpolate(0.0, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new UnseenState().interpolatedColor(20, 100),
        new ColorUtils().interpolate(0.2, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new UnseenState().interpolatedColor(40, 100),
        new ColorUtils().interpolate(0.4, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));

    t.checkExpect(new StartState().interpolatedColor(40, 100),
        new ColorUtils().interpolate(0.4, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new StartState().interpolatedColor(70, 100),
        new ColorUtils().interpolate(0.7, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new StartState().interpolatedColor(40, 100),
        new ColorUtils().interpolate(0.4, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new FinishState().interpolatedColor(40, 100),
        new ColorUtils().interpolate(0.4, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new FinishState().interpolatedColor(70, 100),
        new ColorUtils().interpolate(0.7, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));
    t.checkExpect(new FinishState().interpolatedColor(40, 100),
        new ColorUtils().interpolate(0.4, Constants.CELL_CLOSE_COLOR, Constants.CELL_FAR_COLOR));

    // Only the unseen and finish and start states appear any different
    t.checkExpect(new SeenState().interpolatedColor(10, 100), Constants.CELL_SEEN_COLOR);
    t.checkExpect(new PathState().interpolatedColor(10, 100), Constants.CELL_PATH_COLOR);
  }

  // Tests methods associated with adding, removing, and
  // checking for whether or not an edge is considered
  // to be part of the minimum spanning tree. These methods
  // are observations/have side effects and work in tandem
  void testCellEdgeTreeMethods(Tester t) {
    this.initTestData();
    CellEdge edge = new CellEdge(new DummyCell(), new DummyCell());

    // The edge is not in the tree
    t.checkExpect(edge.inTree(), false);

    edge.becomeEdgeInTree();

    // The edge is now in the tree
    t.checkExpect(edge.inTree(), true);

    edge.becomeEdgeInTree();

    // The edge is still in the tree
    t.checkExpect(edge.inTree(), true);

    edge.becomeEdgeInTree();
    edge.becomeEdgeInTree();
    edge.becomeEdgeInTree();

    // The edge is STILL in the tree
    t.checkExpect(edge.inTree(), true);

    edge.removeFromTree();

    // The edge is no longer in the tree
    t.checkExpect(edge.inTree(), false);

    edge.removeFromTree();
    edge.removeFromTree();
    edge.removeFromTree();

    // The edge is STILL no longer in the tree
    t.checkExpect(edge.inTree(), false);

    edge.becomeEdgeInTree();
    edge.removeFromTree();

    // The edge is STILL no longer in the tree
    t.checkExpect(edge.inTree(), false);
  }

  void testCellEdgeCompareTo(Tester t) {
    this.initTestData();

    // Comparison is a partial order on the
    // set of edges (so we check it has the nice properties
    // that such operators have)

    CellEdge edge1 = new CellEdge(new DummyCell(), new DummyCell(), 10);
    CellEdge edge2 = new CellEdge(new DummyCell(), new DummyCell(), 20);
    CellEdge edge3 = new CellEdge(new DummyCell(), new DummyCell(), 30);
    CellEdge edge4 = new CellEdge(new DummyCell(), new DummyCell(), 40);

    t.checkExpect(edge1.compareTo(edge1) == 0, true); // A <= A
    t.checkExpect(edge1.compareTo(edge2) < 0, true); // A <= B
    t.checkExpect(edge2.compareTo(edge1) > 0, true); // but B <= A not true

    // Testing transitivity (A >= B, B >= C -> A >= C)
    t.checkExpect(edge1.compareTo(edge2) < 0, true);
    t.checkExpect(edge2.compareTo(edge3) < 0, true);
    t.checkExpect(edge1.compareTo(edge3) < 0, true);

    t.checkExpect(edge3.compareTo(edge2) > 0, true);
    t.checkExpect(edge2.compareTo(edge1) > 0, true);
    t.checkExpect(edge3.compareTo(edge1) > 0, true);
  }

  void testCellEdgeFirstConnectionThatIsnt(Tester t) {
    this.initTestData();

    ICell cell1 = new Cell();
    ICell cell2 = new Cell();
    ICell cell3 = new Cell();
    ICell cell4 = new Cell();
    ICell inNothing = new Cell();
    CellEdge edge1 = new CellEdge(cell1, cell2, 10);
    CellEdge edge2 = new CellEdge(cell2, cell1, 20);
    CellEdge edge3 = new CellEdge(cell2, cell3, 9);
    CellEdge edge4 = new CellEdge(cell4, cell4, 40);
    CellEdge edge5 = new CellEdge(cell3, cell3, 0);

    t.checkExpect(edge1.firstConnectionThatIsnt(cell1), cell2);
    t.checkExpect(edge1.firstConnectionThatIsnt(cell2), cell1);
    t.checkExpect(edge1.firstConnectionThatIsnt(inNothing), cell1);
    t.checkExpect(edge2.firstConnectionThatIsnt(cell2), cell1);
    t.checkExpect(edge2.firstConnectionThatIsnt(cell1), cell2);
    t.checkExpect(edge1.firstConnectionThatIsnt(inNothing), cell2);
    t.checkExpect(edge3.firstConnectionThatIsnt(cell2), cell3);
    t.checkExpect(edge3.firstConnectionThatIsnt(cell3), cell2);
    t.checkExpect(edge3.firstConnectionThatIsnt(inNothing), cell3);

    t.checkException(
        new IllegalArgumentException("This CellEdge points to the given cell in both directions."),
        edge4, "firstConnectionThatIsnt", cell4);
    t.checkException(
        new IllegalArgumentException("This CellEdge points to the given cell in both directions."),
        edge5, "firstConnectionThatIsnt", cell3);
  }

  void testICellToImage(Tester t) {
    this.initTestData();

    ICell cell = new Cell();
    ICell hexCell = new HexCell();

    t.checkExpect(cell.toImage(100), this.boxCellImageUnseen100x100);
    t.checkExpect(cell.toImage(150), this.boxCellImageUnseen150x150);
    t.checkExpect(hexCell.toImage(100), this.hexCellImageUnseen100x100);

    // Tell the cell it has been seen
    cell.enter(new SeenState());
    hexCell.enter(new SeenState());

    t.checkExpect(cell.toImage(100), this.boxCellImageSeen100x100);
    t.checkExpect(cell.toImage(150), this.boxCellImageSeen150x150);
    t.checkExpect(hexCell.toImage(100), this.hexCellImageSeen100x100);

    // Tell the cell it is the "first" cell
    cell.enter(new StartState());
    hexCell.enter(new StartState());

    t.checkExpect(cell.toImage(100), this.boxCellImageStart100x100);
    t.checkExpect(cell.toImage(150), this.boxCellImageStart150x150);
    t.checkExpect(hexCell.toImage(100), this.hexCellImageStart100x100);

    // Tell the cell it is the "final" cell
    cell.enter(new FinishState());
    hexCell.enter(new FinishState());

    t.checkExpect(cell.toImage(100), this.boxCellImageFinish100x100);
    t.checkExpect(cell.toImage(150), this.boxCellImageFinish150x150);
    t.checkExpect(hexCell.toImage(100), this.hexCellImageFinish100x100);

    // Tell the cell it is in the final path through the maze
    cell.enter(new PathState());
    hexCell.enter(new PathState());

    t.checkExpect(cell.toImage(100), this.boxCellImagePath100x100);
    t.checkExpect(cell.toImage(150), this.boxCellImagePath150x150);
    t.checkExpect(hexCell.toImage(100), this.hexCellImagePath100x100);
  }

  void testICellToImageInterpolated(Tester t) {
    this.initTestData();

    ICell cell = new Cell();
    ICell hexCell = new HexCell();

    t.checkExpect(cell.toImage(100, 10, 100), this.boxCellImageUnseenInterpolated100x100);
    t.checkExpect(cell.toImage(150, 10, 100), this.boxCellImageUnseenInterpolated150x150);
    t.checkExpect(hexCell.toImage(100, 10, 100), this.hexCellImageUnseenInterpolated100x100);

    // Tell the cell it has been seen
    cell.enter(new SeenState());
    hexCell.enter(new SeenState());

    t.checkExpect(cell.toImage(100, 50, 100), this.boxCellImageSeenInterpolated100x100);
    t.checkExpect(cell.toImage(150, 50, 100), this.boxCellImageSeenInterpolated150x150);
    t.checkExpect(hexCell.toImage(100, 10, 100), this.hexCellImageSeenInterpolated100x100);

    // Tell the cell it is the "first" cell
    cell.enter(new StartState());
    hexCell.enter(new StartState());

    t.checkExpect(cell.toImage(100, 50, 100), this.boxCellImageStartInterpolated100x100);
    t.checkExpect(cell.toImage(150, 50, 100), this.boxCellImageStartInterpolated150x150);
    t.checkExpect(hexCell.toImage(100, 10, 100), this.hexCellImageStartInterpolated100x100);

    // Tell the cell it is the "final" cell
    cell.enter(new FinishState());
    hexCell.enter(new FinishState());

    t.checkExpect(cell.toImage(100, 50, 100), this.boxCellImageFinishInterpolated100x100);
    t.checkExpect(cell.toImage(150, 50, 100), this.boxCellImageFinishInterpolated150x150);
    t.checkExpect(hexCell.toImage(100, 10, 100), this.hexCellImageFinishInterpolated100x100);

    // Tell the cell it is in the final path through the maze
    cell.enter(new PathState());
    hexCell.enter(new PathState());

    t.checkExpect(cell.toImage(100, 50, 100), this.boxCellImagePathInterpolated100x100);
    t.checkExpect(cell.toImage(150, 50, 100), this.boxCellImagePathInterpolated150x150);
    t.checkExpect(hexCell.toImage(100, 10, 100), this.hexCellImagePathInterpolated100x100);
  }

  void testOurStack(Tester t) {
    this.initTestData();
    ICollection<String> stackOStrings = new OurStack<String>();
    ICollection<Integer> stackOInts = new OurStack<Integer>();

    // Check there is nothing in the stacks
    t.checkExpect(stackOStrings.isEmpty(), true);
    t.checkExpect(stackOInts.isEmpty(), true);

    // Add some items
    stackOStrings.add("More");
    stackOStrings.add("Fun");
    stackOStrings.add("In");
    stackOStrings.add("Fundies");

    stackOInts.add(10);
    stackOInts.add(-100);
    stackOInts.add(0);
    stackOInts.add(9);

    // There should be items in the stack now
    t.checkExpect(stackOStrings.isEmpty(), false);
    t.checkExpect(stackOInts.isEmpty(), false);

    // Now we can check that removing works as we expect
    t.checkExpect(stackOStrings.remove(), "Fundies");
    t.checkExpect(stackOStrings.remove(), "In");
    t.checkExpect(stackOStrings.remove(), "Fun");
    t.checkExpect(stackOStrings.remove(), "More");

    t.checkExpect(stackOInts.remove(), 9);
    t.checkExpect(stackOInts.remove(), 0);
    t.checkExpect(stackOInts.remove(), -100);
    t.checkExpect(stackOInts.remove(), 10);

    // Attempting to remove more items when there are none
    // is an exception
    t.checkException(new EmptyStackException(), stackOStrings, "remove");
    t.checkException(new EmptyStackException(), stackOInts, "remove");
  }

  void testOurQueue(Tester t) {
    this.initTestData();
    ICollection<String> queueOStrings = new Queue<String>();
    ICollection<Integer> queueOInts = new Queue<Integer>();

    // Check there is nothing in the queues
    t.checkExpect(queueOStrings.isEmpty(), true);
    t.checkExpect(queueOInts.isEmpty(), true);

    // Add some items
    queueOStrings.add("More");
    queueOStrings.add("Fun");
    queueOStrings.add("In");
    queueOStrings.add("Fundies");

    queueOInts.add(10);
    queueOInts.add(-100);
    queueOInts.add(0);
    queueOInts.add(9);

    // There should be items in the stack now
    t.checkExpect(queueOStrings.isEmpty(), false);
    t.checkExpect(queueOInts.isEmpty(), false);

    // Now we can check that removing works as we expect
    t.checkExpect(queueOStrings.remove(), "More");
    t.checkExpect(queueOStrings.remove(), "Fun");
    t.checkExpect(queueOStrings.remove(), "In");
    t.checkExpect(queueOStrings.remove(), "Fundies");

    t.checkExpect(queueOInts.remove(), 10);
    t.checkExpect(queueOInts.remove(), -100);
    t.checkExpect(queueOInts.remove(), 0);
    t.checkExpect(queueOInts.remove(), 9);

    // Attempting to remove more items when there are none
    // is an exception
    t.checkException(new NoSuchElementException(), queueOStrings, "remove");
    t.checkException(new NoSuchElementException(), queueOInts, "remove");
  }

  void testNotPred(Tester t) {
    this.initTestData();
    t.checkExpect(new NotPred<Posn>(new InFirstQuadrant()).test(new Posn(-100, 10)), true);
    t.checkExpect(new NotPred<Integer>(new IsPositive()).test(10), false);
    t.checkExpect(new NotPred<Integer>(new IsPositive()).test(-10), true);
    t.checkExpect(new NotPred<Integer>(new IsPositive()).test(-129), true);
  }

  void testMazeWorldOnTick(Tester t) {
    this.initTestData();

    // this maze world should have a maze identical to the starting maze
    MazeWorld w = new MazeWorld(3, 3, 90, 90, 0, false);
    w.onKeyEvent("!");

    w.onTick();

    // the starting world after the maze is drawn, where we ared at the start
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // take a step
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // one more step
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // see nothing change when we pause
    w.onKeyEvent(" ");
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));
    w.onKeyEvent(" ");

    // a few more steps, so we can see the thing when it is finished
    w.onTick();
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_PATH_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_PATH_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_PATH_COLOR, Constants.CELL_PATH_COLOR,
            Constants.CELL_PATH_COLOR));

    // now we'll start a breadth first search
    w.onKeyEvent("b");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // this does not encompass all possible onKey options, but it does go
    // through each step of the onTick process. Further usage of onTick with
    // the remaining onKey options will be used in the onKeyEvent tests.
  }

  void testMazeWorldOnKeyEvent(Tester t) {
    this.initTestData();

    // begin with the example 3x3 maze

    MazeWorld w = new MazeWorld(3, 3, 90, 90, 0, false);

    // EVENT 1: break all walls
    w.onKeyEvent("!");

    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // finish this maze
    w.onTick();
    w.onTick();
    w.onTick();
    w.onTick();

    // EVENT 2: breadth first search
    w.onKeyEvent("b");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // make one more step, to make sure it's breadth first

    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FUTUREHEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // EVENT 3: hide all seen cells

    w.onKeyEvent("-");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    w.onKeyEvent("-");

    // EVENT 3: depth first search

    w.onKeyEvent("p");
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 4: pause

    w.onKeyEvent(" ");
    w.onTick();
    // make sure nothing happens
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    w.onKeyEvent(" ");

    // take another step, to make sure it's depth first

    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 5: manual search
    w.onKeyEvent("m");
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // make sure that it does not run automatically
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 6: move down on manual
    w.onKeyEvent("down");

    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_HEAD_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 7: move right on manual

    w.onKeyEvent("right");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_HEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    w.onKeyEvent("right");
    w.onTick();

    // EVENT 8: move up on manual
    w.onKeyEvent("up");
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_HEAD_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
            Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 9: move left on manual
    w.onKeyEvent("left");
    w.onTick();
    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors(Constants.CELL_SEEN_COLOR, Constants.CELL_HEAD_COLOR,
            Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR, Constants.CELL_SEEN_COLOR,
            Constants.CELL_SEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 10: make a new unbiased maze
    w.onKeyEvent("n");
    w.onTick();
    w.onKeyEvent("!");

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColors2(Constants.CELL_START_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_FINISH_COLOR));

    // EVENT 11: make a new rowed maze
    w.onKeyEvent("r");
    w.onTick();
    w.onKeyEvent("!");

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColorsRow(Constants.CELL_START_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // EVENT 12: show all distances from the start
    w.onKeyEvent("k");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColorsRow(Constants.CELL_HEAD_COLOR,
            this.cellColorInterpolated(1, 8), this.cellColorInterpolated(2, 8),
            this.cellColorInterpolated(5, 8), this.cellColorInterpolated(4, 8),
            this.cellColorInterpolated(3, 8), this.cellColorInterpolated(6, 8),
            this.cellColorInterpolated(7, 8), this.cellColorInterpolated(8, 8)));

    w.onKeyEvent("k");

    // EVENT 13: show all distances from the finish
    w.onKeyEvent("l");
    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColorsRow(Constants.CELL_HEAD_COLOR,
            this.cellColorInterpolated(7, 8), this.cellColorInterpolated(6, 8),
            this.cellColorInterpolated(3, 8), this.cellColorInterpolated(4, 8),
            this.cellColorInterpolated(5, 8), this.cellColorInterpolated(2, 8),
            this.cellColorInterpolated(1, 8), this.cellColorInterpolated(0, 8)));

    w.onKeyEvent("l");

    // EVENT 14: make a new columned maze
    w.onKeyEvent("c");
    w.onTick();
    w.onKeyEvent("!");

    t.checkExpect(w.makeScene(),
        this.drawExample3x3MazeWithColorsCol(Constants.CELL_START_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // EVENT 15: Make a hexagonal maze (a new world to test this)

    w.onKeyEvent("q");
    w.onTick();
    w.onKeyEvent("!");

    t.checkExpect(w.makeScene(),
        this.drawExample3x3HexagonMazeWithColors2(Constants.CELL_START_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    // New world for testing hexagons further
    w = new MazeWorld(3, 3, 800, 800, 0, true);

    w.onKeyEvent("!");

    t.checkExpect(w.makeScene(),
        this.drawExample3x3HexagonMazeWithColors(Constants.CELL_START_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));

    w.onTick();

    t.checkExpect(w.makeScene(),
        this.drawExample3x3HexagonMazeWithColors(Constants.CELL_HEAD_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR, Constants.CELL_UNSEEN_COLOR,
            Constants.CELL_UNSEEN_COLOR, Constants.CELL_FINISH_COLOR));
  }

  void testWouldContainInBox(Tester t) {
    this.initTestData();

    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(5.0, new Point2D(12.0, 10.0)), true);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(5.0, new Point2D(14.0, 10.0)), true);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(5.0, new Point2D(6.0, 9.0)), true);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(5.0, new Point2D(4.0, 10.0)), false);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(5.0, new Point2D(9.0, 10.0)), true);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(1.0, new Point2D(12.0, 10.0)), false);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(7.0, new Point2D(0.0, 10.0)), false);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(9.0, new Point2D(12.0, 10.0)), true);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(1.5, new Point2D(12.0, 10.0)), false);
    t.checkExpect(new Point2D(10.0, 10.0).wouldContainInBox(2.5, new Point2D(12.0, 10.0)), true);
    t.checkExpect(new Point2D(0.0, 3.0).wouldContainInBox(10.5, new Point2D(20.0, 12.0)), false);
  }

  void testConvertToWorld(Tester t) {
    this.initTestData();

    t.checkExpect(new StandardCoordinateSystem().convertToWorld(new Point2D(0.0, 3.0))
        .wouldContainInBox(0.01, new Point2D(0.0, 3.0)), true);
    t.checkExpect(new StandardCoordinateSystem().convertToWorld(new Point2D(1.0, 9.0))
        .wouldContainInBox(0.01, new Point2D(1.0, 9.0)), true);
    t.checkExpect(new StandardCoordinateSystem().convertToWorld(new Point2D(0.0, 3.0))
        .wouldContainInBox(0.01, new Point2D(0.0, 3.0)), true);
    t.checkExpect(new GraphicsCoordinateSystem().convertToWorld(new Point2D(12.5, 3.0))
        .wouldContainInBox(0.01, new Point2D(12.5, -3.0)), true);
    t.checkExpect(new GraphicsCoordinateSystem().convertToWorld(new Point2D(-100.4, -30.4))
        .wouldContainInBox(0.01, new Point2D(-100.4, 30.4)), true);
    t.checkExpect(
        new RotatedCoordinateSystem(0.0, new StandardCoordinateSystem())
            .convertToWorld(new Point2D(12.5, 3.0)).wouldContainInBox(0.01, new Point2D(12.5, 3.0)),
        true);

    t.checkExpect(new RotatedCoordinateSystem(Math.PI / 2, new StandardCoordinateSystem())
        .convertToWorld(new Point2D(80.0, -190.0))
        .wouldContainInBox(0.01, new Point2D(190.0, 80.0)), true);
    t.checkExpect(new RotatedCoordinateSystem(Math.PI, new StandardCoordinateSystem())
        .convertToWorld(new Point2D(80.0, -190.0))
        .wouldContainInBox(0.01, new Point2D(-80.0, 190.0)), true);
  }

  void testConvertToThis(Tester t) {
    this.initTestData();

    t.checkExpect(new StandardCoordinateSystem().convertToThis(new Point2D(0.0, 3.0))
        .wouldContainInBox(0.01, new Point2D(0.0, 3.0)), true);
    t.checkExpect(new StandardCoordinateSystem().convertToThis(new Point2D(1.0, 9.0))
        .wouldContainInBox(0.01, new Point2D(1.0, 9.0)), true);
    t.checkExpect(new StandardCoordinateSystem().convertToThis(new Point2D(0.0, 3.0))
        .wouldContainInBox(0.01, new Point2D(0.0, 3.0)), true);
    t.checkExpect(new GraphicsCoordinateSystem().convertToThis(new Point2D(12.5, 3.0))
        .wouldContainInBox(0.01, new Point2D(12.5, -3.0)), true);
    t.checkExpect(new GraphicsCoordinateSystem().convertToThis(new Point2D(-100.4, -30.4))
        .wouldContainInBox(0.01, new Point2D(-100.4, 30.4)), true);
    t.checkExpect(
        new RotatedCoordinateSystem(0.0, new StandardCoordinateSystem())
            .convertToThis(new Point2D(12.5, 3.0)).wouldContainInBox(0.01, new Point2D(12.5, 3.0)),
        true);

    t.checkExpect(new RotatedCoordinateSystem(Math.PI / 2, new StandardCoordinateSystem())
        .convertToThis(new Point2D(80.0, -190.0))
        .wouldContainInBox(0.01, new Point2D(-190.0, -80.0)), true);
    t.checkExpect(new RotatedCoordinateSystem(Math.PI, new StandardCoordinateSystem())
        .convertToThis(new Point2D(80.0, -190.0))
        .wouldContainInBox(0.01, new Point2D(-80.0, 190.0)), true);
  }
}