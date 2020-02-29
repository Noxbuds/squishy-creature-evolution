import java.util.Collections;
import java.util.Comparator;
import java.util.Arrays;

/////////////////////////////////////////////
//                                         //
//               CONSTANTS                 //
//                                         //
/////////////////////////////////////////////

// Window dimensions
final int windWidth = 1000;
final int windHeight = 1000;

// World to screen multiplier
final float worldToScreenMult = 100f;

// Viewport dimensions
final float viewWidth = windWidth / worldToScreenMult;
final float viewHeight = windHeight / worldToScreenMult;

// RNG seed
final int seed = 1234;

// Acceleration due to gravity
final float gravity = 9.81;

// Aerodynamic drag
final float fluidDensity = 1.225;

// Extension amount for springs
final float extendAmount = 1.5;

// Extension threshold
final float extendThreshold = 0.01;

// Extension timestep
final float extendTimestep = 0.15;

// Creature squishiness
final float creatureSquishiness = 1500.0;

// Whether to save
final boolean saveGenerations = false;

// Simulation timestep
final float simTimestep = 1.0 / 300.0;

// Simulation time (seconds, usually 15)
final float simTime = 15.0;

// Maximum number of generations to auto-sim
// When the app starts, it will automatically begin simulating
// until it reaches this number of generations. However, if you
// set it to 0, it will not automatically simulate and you will
// be able to manually control it
final int maxAutoGens = 0;

// Whether muscles instantly contract/expand or gradually change
final boolean instantMuscleExpand = false;

//
//
// Creature start properties
//
//

// Starting size
final float startCellSize = 0.5;

// Grid width and height
final int startGridWidth = 4;
final int startGridHeight = 4;

// Start position
final PVector startPosition = new PVector(-(startGridWidth + 1) * startCellSize * 2.0, 0.2);

// Number of creatures
final int creatureCount = 500;

// Mutation chance
final float mutateChance = 0.01;

// Pre-simulate time
final float presimTime = 0.1;

// Node mass
final float nodeMass = 1;

/////////////////////////////////////////////
//                                         //
//            GLOBAL VARIABLES             //
//                                         //
/////////////////////////////////////////////

// Creature view type
//
// 0 - normal (cells, no nodes)
// 1 - cell + node
// 2 - wireframe
//
int creatureViewType = 0;

// Whether we're currently simulating
boolean simulating;
boolean batchInProgress;

// Generation count
int generationCount;

// Pre-sim timer for waiting
float presimTimer;

// View offset
float viewOffset;

// Creatures saved at certain percentile
SquishyCreature creature0th;
SquishyCreature creature50th;
SquishyCreature creature100th;

/////////////////////////////////////////////
//                                         //
//               UTILITIES                 //
//                                         //
/////////////////////////////////////////////

// Maps a world x co-ordinate to screen space
float mapX(float x)
{
  return x / viewWidth * windWidth;
}

// Maps a world y co-ordinate to screen space
float mapY(float y)
{
  return (viewHeight - y) / viewHeight * windHeight;
}

// Maps screen co-ordinates to world co-ordinates
float screenToWorldX(float pos)
{
  return (pos / windWidth) * viewWidth;
}

// Converts a cell type to an ID
int cellTypeToID(char type)
{
  switch (type)
  {
    default:
    case '_':
      return 0;
    case 'E':
      return 1;
    case 'C':
      return 2;
    case '#':
      return 3;
  }
}

// Converts a cell ID to a type
char cellIDToType(int id)
{
  switch (id)
  {
    default:
    case 0:
      return '_';
    case 1:
      return 'E';
    case 2:
      return 'C';
    case 3:
      return '#';
  }
}

/////////////////////////////////////////////
//                                         //
//                CLASSES                  //
//                                         //
/////////////////////////////////////////////

// Stores data about a node
//
// - Each node has a radius; this should likely just be fixed
// - Each node will have a position and velocity
// - Each node will have a list of forces that it will resolve on the next frame
//
class BouncyNode
{
  // Radius and mass
  public float radius;
  public float mass;
  
  // Position
  public PVector position;
  
  // Velocity
  public PVector velocity;
  
  // Instantaneous acceleration
  private PVector instAccel;
  
  // Drag coefficient
  private float dragCoeff;
  
  // Friction coefficient
  private float frictionCoeff;
  
  // Whether the node is on the outside of a creature
  private boolean outsideNode;
  
  // Forces to resolve next frame
  private ArrayList<PVector> forces;
  
  // Constructor for a bouncy node
  public BouncyNode(float radius, PVector position, float dragCoeff, float frictionCoeff, boolean outsideNode)
  {
    // Set properties
    this.radius = radius;
    this.mass = nodeMass;
    this.position = position;
    this.dragCoeff = dragCoeff;
    this.frictionCoeff = frictionCoeff;
    this.outsideNode = outsideNode;
    
    // Set velocity to zero
    this.velocity = new PVector(0, 0);
    
    // Initialise forces array
    forces = new ArrayList<PVector>();
  }
  
  // Adds a force to the list
  public void AddForce(PVector force)
  {
    forces.add(force);
  }
  
  // Calculates air resistance
  private void CalculateAirResistance()
  {
    // Air resistance
    // Fd = 1/2 p v^2 Cd A
    //
    // - p is the density of the fluid (air = 1.225)
    // - v is the object's speed
    // - A is the cross-sectional area (2d so area is half circumference)
    // - Cd is drag coefficient (0.5 for a circle?)
    float p = fluidDensity;
    float v2 = velocity.magSq();
    float A = 6.28318 * radius; // (2pi times radius)
    
    // Get a direction for the air resistance vector
    PVector u = PVector.mult(velocity, -1.0);
    u.normalize();
    
    // Air resistance
    PVector Fd = PVector.mult(u, 0.5 * p * v2 * dragCoeff * A);
    
    // Add the air resistance
    AddForce(Fd);
  }
  
  // Calculates the normal contact force
  private void CalculateContactForce()
  { 
    // Virtual extension/squishing of the node
    float extension = position.y - radius;
    
    // Spring extension constant
    float k = 2000;
    
    // Setup contact force
    PVector N = new PVector(0, 0);
    
    // Check if we're colliding with the ground
    if (extension < 0)
    {
      // Contact force, velocity resisting force, and spring force
      N.y = gravity * mass;
      N.y += k * 0.1 * -velocity.y;
      N.y += (radius * 0.95 - position.y) * k;
      
      // Setup friction
      // F = -uRv
      //
      // - u is friction coefficient
      // - R is magnitude of contact force
      // - v is unit vector of velocity
      PVector v = velocity.copy();
      v.normalize();
      
      PVector friction = PVector.mult(v, -frictionCoeff * abs(N.y));
      friction.y = 0; // just x-axis friction
      
      // Apply friction and normal force
      AddForce(N);
      AddForce(friction);
    }
  }
  
  // Resolves all forces
  public void ResolveForces(float timestep)
  {
    // Reset instantaneous acceleration
    instAccel = new PVector(0, 0);
    
    // Loop over each force
    for (int i = forces.size() - 1; i >= 0; i--)
    {
      // Calculate acceleration
      // F = ma  ==>  a = F/m
      PVector acceleration = PVector.div(forces.get(i), mass);
      
      // Add the acceleration
      instAccel.add(acceleration);
      
      // Add the acceleration, multiplied by time step, to the velocity
      velocity.add(PVector.mult(acceleration, timestep));
      
      // Remove this force from the list
      forces.remove(i);
    }
  }
  
  // Calculates a physics timestep
  // External forces should be added first
  public void PhysicsStep(float timestep)
  {
    // Calculate the contact force
    CalculateContactForce();
    
    // Do air resistance
    if (outsideNode)
      CalculateAirResistance();
    
    // Resolve forces
    ResolveForces(timestep);
    
    // Move the position
    position.add(PVector.mult(velocity, timestep));
  }
}

// Contains a pair of nodes which are connected via a spring
class SpringPair
{
  // The nodes in the pair
  private BouncyNode nodeA;
  private BouncyNode nodeB;
  
  // The spring constant
  private float springConstant;
  
  // The ideal length
  private float initialLength;
  private float targetLength;
  
  // Extend/contract stuff
  private float extendTimer;
  private float extendInterval;
  private boolean extended;
  
  // Whether it is actually a muscle
  private boolean isMuscle;
  
  // Constructor
  public SpringPair(BouncyNode nodeA, BouncyNode nodeB, float springConstant, float targetLength, boolean startExtended, int extendInterval, boolean isMuscle)
  {
    this.nodeA = nodeA;
    this.nodeB = nodeB;
    this.springConstant = springConstant;
    this.initialLength = targetLength;
    this.targetLength = startExtended ? targetLength * extendAmount : targetLength;
    this.extended = extendInterval < 0.1 ? false : startExtended;
    this.extendInterval = extendInterval;
    this.isMuscle = isMuscle;
  }
  
  // Calculates the forces on the pair
  public void PhysicsStep(float timestep)
  {
    // Increment extend timer
    extendTimer += timestep;
    
    // Check if this is a muscle (only muscles expand/contract)
    if (isMuscle)
    {
      // Check if the extend timer exceeds the interval
      if (extendTimer >= extendInterval * extendTimestep && extendInterval > 0)
      {
        // Reset the timer
        extendTimer = 0;
        
        // Check if the spring is extended or not
        if (extended) extended = false;
        else extended = true;
      }
      
      // Check what type of contraction/extension we are after
      if (instantMuscleExpand)
      {
        if (extended && extendInterval > 0) targetLength = initialLength * extendAmount;
        else targetLength = initialLength;
      }
      else
      {
        // Mix the target length
        if (extended) targetLength = lerp(initialLength, initialLength * extendAmount, extendTimer / (extendInterval * extendTimestep));
        else if (extendInterval > 0) targetLength = lerp(initialLength * extendAmount, initialLength, extendTimer / (extendInterval * extendTimestep));
      }
    }
    
    // Vector pointing from A to B
    PVector u = PVector.sub(nodeA.position, nodeB.position);
    
    // Calculate length between the nodes
    float len = u.mag();
    
    // Calculate extension
    float extension = len - targetLength;
    
    // If extension is below a threshold, ignore it
    if (abs(extension) < extendThreshold) return;
    
    // Normalize the direction vector
    u.normalize();
    
    // F = -ukx
    // - u is the directional vector
    // - k is the spring constant
    // - x is the extension
    PVector force = PVector.mult(u, springConstant * extension);
    
    // Apply the force to each node
    nodeA.AddForce(PVector.mult(force, -1.0));
    nodeB.AddForce(force);
  }
}

// This class stores data about a creature made from squishy cells
// that extend or contract on individual timers. Making use of springs
// to give the illusion of a squishy creature makes them much more believable.
// Each creature will be made up of a grid of cells, and each corner of a cell is a node.
class SquishyCreature implements Comparable<SquishyCreature>
{ 
  // Fitness
  public float fitness;
  
  // The world position of the top-left corner
  private PVector position;
  
  // Grid size
  private int gridWidth;
  private int gridHeight;
  
  // Unit length
  private float cellSize;
  private float cellDiagonalSize;
  
  // Squishiness
  private float squishiness;
  
  // The creature's DNA
  private String DNA;
  
  // All the nodes of the square
  private BouncyNode[] nodes;
  
  // All the node pairs
  private ArrayList<SpringPair> pairs;
  
  // Constructor. Takes in a string for DNA, cell size, and grid size.
  //
  // DNA formatting example (3x3 grid):
  // "EC#_##_EEE25____623"
  // - 3x3 grid unwrapped into 2d, first 3 = first row, second 3 = second row etc
  // - last 6 values are the amount of time steps it takes to extend/contract
  // - E makes it extend first
  // - C makes it contract first
  // - # is just tissue (still squishy though)
  // - _ is an empty space
  public SquishyCreature(PVector position, String DNA, float cellSize, int gridWidth, int gridHeight)
  {
    // Set the base properties
    this.position = position;
    this.DNA = DNA;
    this.cellSize = cellSize;
    this.cellDiagonalSize = 1.41421356 * cellSize;
    this.gridWidth = gridWidth;
    this.gridHeight = gridHeight;
    this.squishiness = creatureSquishiness;
    
    // Initialise lists
    nodes = new BouncyNode[(gridWidth + 1) * (gridHeight + 1)];
    pairs = new ArrayList<SpringPair>();
    
    // Calculate the offset for the frame counters
    int frameCountOffset = gridWidth * gridHeight;
    
    // Loop over each node in the grid
    for (int y = 0; y < gridHeight + 1; y++)
    {
      for (int x = 0; x < gridWidth + 1; x++)
      {
        // Add a node at this point
        nodes[x + y * (gridWidth + 1)] = new BouncyNode(0.05, PVector.add(position, new PVector(x * cellSize, y * cellSize)), 2.0, 1.0, true);
      } 
    }
    
    // Loop over each cell
    for (int y = 0; y < gridHeight; y++)
    {
      for (int x = 0; x < gridWidth; x++)
      {
        // Cell ID
        int cellID = x + y * gridWidth;
        
        // Get cell type
        char cellType = DNA.charAt(cellID);
        
        // Skip if cell type is empty
        if (cellType == '_') continue;
        
        // Extended interval
        int extendInterval = Character.getNumericValue(DNA.charAt(cellID + frameCountOffset));
        
        // If it's a piece of solid tissue, make sure it doesn't extend
        if (cellType == '#') extendInterval = -1;
        
        // If the extend interval is zero, set it to -1 just in case
        if (extendInterval == 0) extendInterval = -1;
        
        // Add the square
        AddSquare(x, y, cellType, extendInterval);
      }
    }
  }
  
  // Sets a spring pair up based on two IDs
  private void CreatePair(int id1, int id2, float len, boolean startExtended, int extendInterval, float toughness, boolean isMuscle)
  {
    // Add each pair
    pairs.add(new SpringPair(nodes[id1], nodes[id2], toughness, len, startExtended, extendInterval, isMuscle));
  }
  
  // Adds a cell at a position
  private void AddSquare(int x, int y, char cellType, int extendInterval)
  {
    // Start extended?
    boolean startExtended = cellType == 'E' ? true : false;
    
    // Toughness
    float toughness = cellType == '#' ? squishiness * 2.0 : squishiness;
    
    // Get the corner IDs
    int BL = x + y * (gridWidth + 1);
    int BR = (x + 1) + y * (gridWidth + 1);
    int TL = x + (y + 1) * (gridWidth + 1);
    int TR = (x + 1) + (y + 1) * (gridWidth + 1);
    
    // We are going from bottom-left to top-right when building the model
    // For most cells, we just want the left and bottom nodes linked as passive links
    CreatePair(BL, BR, cellSize, startExtended, extendInterval, toughness, false);
    CreatePair(BL, TL, cellSize, startExtended, extendInterval, toughness, false);
    
    // Cell IDs used for checking cells to the right and top
    int rightCellID = x + 1 + y * gridWidth;
    int topCellID = x + (y + 1) * gridWidth;
    
    // However, for the right-most and top-most cells, we also want to link their respective edges
    // We will check if we are either at the edge of the grid, or if the right or top cell is empty
    boolean rightCheck = x == gridWidth - 1;
    if (x < gridWidth - 1)
      rightCheck = DNA.charAt(rightCellID) == '_';
    
    boolean topCheck = y == gridHeight - 1;
    if (y < gridHeight - 1)
      topCheck = DNA.charAt(topCellID) == '_';
    
    // Create the pairs mentioned above as required
    if (rightCheck) CreatePair(BR, TR, cellSize, startExtended, extendInterval, toughness, false);
    if (topCheck) CreatePair(TL, TR, cellSize, startExtended, extendInterval, toughness, false);
    
    // Create the diagonal muscle links
    CreatePair(BL, TR, cellDiagonalSize, startExtended, extendInterval, toughness * 2, true);
    CreatePair(TL, BR, cellDiagonalSize, startExtended, extendInterval, toughness * 2, true);
  }
  
  // Does a physics timestep
  public void PhysicsStep(float timestep)
  {
    // Loop over each spring pair
    for (int i = 0; i < pairs.size(); i++)
      pairs.get(i).PhysicsStep(timestep);
    
    // Loop over each node
    for (int i = 0; i < nodes.length; i++)
    {
      // Make sure node is not null
      if (nodes[i] != null)
      {
        // Add gravity
        nodes[i].AddForce(new PVector(0, -gravity));
        
        // Calculate physics
        nodes[i].PhysicsStep(timestep);
        
        // Check if it is underground
        if (nodes[i].position.y < -2.0)
          nodes[i] = null;
      }
    }
  }
  
  // Draws the creature on-screen
  public void Draw()
  {
    // Draw the cells
    if (creatureViewType >= 0)
    {
      // Loop over the grid
      for (int y = 0; y < gridHeight; y++)
      {
        for (int x = 0; x < gridWidth; x++)
        {
          // Extended grid width
          int extendedGridWidth = gridWidth + 1;
          
          // Calculate the IDs for each corner
          int idBL = x + y * extendedGridWidth;
          int idBR = (x + 1) + y * extendedGridWidth;
          int idTL = x + (y + 1) * extendedGridWidth;
          int idTR = (x + 1) + (y + 1) * extendedGridWidth;
          
          // Get the cell type at this point
          char cellType = DNA.charAt(x + y * gridWidth);
          
          // Stop if empty cell
          if (cellType == '_') continue;
          
          // Cell color
          int cellR = 0;
          int cellG = 0;
          int cellB = 0;
          
          // Set color based on cell type
          if (cellType == '#') { cellR += 100; cellG += 100; cellB += 100; };
          if (cellType == 'E') { cellR += 100; cellG += 200; cellB += 100; };// fill(100 + chargeColor, 200 + chargeColor, 100 + chargeColor); // green
          if (cellType == 'C') { cellR += 200; cellG += 100; cellB += 100; };//fill(200 + chargeColor, 200 + chargeColor, 200 + chargeColor); // red
          
          // Set fill and stroke to cell color
          fill(cellR, cellG, cellB);
          stroke(cellR, cellG, cellB);
          
          // Begin the shape
          beginShape();
          
          // Fetch each corner
          BouncyNode BL = nodes[idBL];
          BouncyNode BR = nodes[idBR];
          BouncyNode TL = nodes[idTL];
          BouncyNode TR = nodes[idTR];
          
          // Add vertices in clockwise order
          if (TL != null) vertex(mapX(TL.position.x), mapY(TL.position.y));
          if (TR != null) vertex(mapX(TR.position.x), mapY(TR.position.y));
          if (BR != null) vertex(mapX(BR.position.x), mapY(BR.position.y));
          if (BL != null) vertex(mapX(BL.position.x), mapY(BL.position.y));
          if (TL != null) vertex(mapX(TL.position.x), mapY(TL.position.y));
          
          // End the shape
          endShape();
        }
      }
    }
    
    // Draw the nodes
    if (creatureViewType >= 1)
    {
      // Loop over each node
      for (int i = 0; i < nodes.length; i++)
      {
        BouncyNode node = nodes[i];
        ellipse(mapX(node.position.x), mapY(node.position.y), node.radius * 2.0 * worldToScreenMult, node.radius * 2.0 * worldToScreenMult);
      }
    }
    
    // Draw wireframe
    if (creatureViewType >= 2)
    {
      // Loop over each pair
      for (int i = 0; i < pairs.size(); i++)
      {
        BouncyNode nodeA = pairs.get(i).nodeA;
        BouncyNode nodeB = pairs.get(i).nodeB;
        line(mapX(nodeA.position.x), mapY(nodeA.position.y), mapX(nodeB.position.x), mapY(nodeB.position.y));
      }
    }
  }
  
  // Copies a creature and mutates it
  public SquishyCreature Copy()
  {
    // New DNA
    String newDNA = "";
    
    // Calculate the offset for the timers
    int timerOffset = gridWidth * gridHeight;
    
    // Loop over letters
    for (int i = 0; i < timerOffset; i++)
    {
      // Roll the dice
      float num = random(0, 1);
      
      // Check if we're gonna mutate
      if (num < mutateChance)
      {
        // Convert cell type to an ID
        int id = cellTypeToID(DNA.charAt(i));
        
        // Pick whether we'll be mutating up or down
        int mutateUp = int(random(0, 1));
        
        // Mutate the ID
        id += mutateUp == 1 ? 1 : -1;
        
        // Convert cell Id to type
        char cellType = cellIDToType(id);
        
        // Check if the new cell is C/E, and if it is make sure
        // that it doesn't have a timer of 0 (helps prevent twitch
        // sliding)
        if (cellType == 'C' || cellType == 'E')
          if (Character.getNumericValue(DNA.charAt(i + timerOffset)) < 1)
            DNA = DNA.substring(0, i + timerOffset) + "1" + DNA.substring(i + timerOffset + 1);
        
        // Add this to the string
        newDNA = newDNA + cellType;
      }
      else
        newDNA = newDNA + DNA.charAt(i);
    }
    
    // Loop over numbers
    for (int i = timerOffset; i < DNA.length(); i++)
    {
      // Roll the dice
      float num = random(0, 1);
      
      // Check if we're gonna mutate
      if (num < mutateChance )
      {
        // Get the number as an integer
        int value = Character.getNumericValue(DNA.charAt(i));
        
        // Pick whether we'll mutate up or down
        int mutateUp = int(random(0, 1));
        
        // Mutate the number
        value += mutateUp == 1 ? 1 : -1;
        
        // Make sure we don't reach zero or change cell type
        if (value < 0) value = 0;
        if (value < 1 && DNA.charAt(i - timerOffset) != '#') value = 1;
        
        // Also make sure we don't exceed 9
        if (value > 9) value = 9;
        
        // Add to the string
        newDNA = newDNA + value;
      }
      else
        newDNA = newDNA + DNA.charAt(i);
    }
    
    // Return the new DNA
    return CreatureFromDNA(newDNA);
  }
  
  // In order to compare creatures, we want to be able to use
  // Java comparison operators (<>). Overriding compareTo(..) will
  // allow us to do this
  public int compareTo(SquishyCreature creature)
  {
    // Just compare fitness
    return Float.compare(this.fitness, creature.fitness);
  }
}

// Randomly generates a creature's DNA
String GenerateNewDNA(int gridWidth, int gridHeight)
{
  // First and second parts of the string
  String firstPart = "";
  String secondPart = "";
  
  // Loop over the grid
  for (int i = 0; i < gridWidth * gridHeight; i++)
  {
    // Pick a number from 0-3 to represent the cell type
    int typeInt = int(random(0, 4));
    
    // Convert the integer to a character
    // 0 = _
    // 1 = E
    // 2 = C
    // 3 = #
    String typeChar = typeInt == 0 ? "_" : (typeInt == 1 ? "E" : (typeInt == 2 ? "C" : "#"));
    
    // Append the character to the end of the string
    firstPart = firstPart + typeChar;
    
    // Pick a random number from 1 to 9 for the extend timer
    int extendInterval = typeChar == "#" ? 0 : int(random(1, 10));
    
    // Append the extend timer to the second part
    secondPart = secondPart + "" + extendInterval;
  }
  
  // Return the second part appended to the first part
  return firstPart + secondPart;
}

// Generates a new creature
SquishyCreature GenerateNewCreature()
{
  // Generate the DNA
  String DNA = GenerateNewDNA(startGridWidth, startGridHeight);
  
  // Create the creature
  return new SquishyCreature(startPosition, DNA, startCellSize, startGridWidth, startGridHeight);
}

// Generates a creature from DNA
SquishyCreature CreatureFromDNA(String DNA)
{
  return new SquishyCreature(startPosition, DNA, startCellSize, startGridWidth, startGridHeight);
}

/////////////////////////////////////////////
//                                         //
//             USER INTERFACE              //
//                                         //
/////////////////////////////////////////////

// Draws the general user interface
void drawUI()
{
  // Set the background to a light grey
  background(220);
  
  // Check the 100% creature exists
  if (creature100th != null)
  {
    // Show the 100% creature in top left
    scale(0.333333);
    fill(220);
    rect(0, 0, 3000, 1000);
    fill(0);
    
    // Generation title
    textSize(110);
    text("Generation " + generationCount, 20, 130);
    
    // Creature
    textSize(60);
    text("1. [[" + creature100th.DNA + "]]   " + creature100th.fitness + " metres", 950, 80);
    //translate((viewWidth / 2 - creature100th.fitness) * worldToScreenMult, -250.0);
    translate(600, 0.0);
    creature100th.Draw();
    translate(-600, 0.0);
    //translate((viewWidth / 2 - creature100th.fitness) * worldToScreenMult * -1, 250.0);
    scale(3.0);
  }
  
  // Check the 50% creature exists
  if (creature50th != null)
  {
    // Show the 0% creature in top left
    translate(0, 333);
    scale(0.3333333);
    fill(180);
    rect(0, 0, 3000, 1000);
    fill(0);
    textSize(60);
    text((creatureCount / 2) + ". [[" + creature50th.DNA + "]]   " + creature50th.fitness + " metres", 950, 80);
    //translate((viewWidth / 2 - creature50th.fitness) * worldToScreenMult, -250.0);
    translate(600, 0.0);
    creature50th.Draw();
    translate(-600, 0.0);
    //translate((viewWidth / 2 - creature50th.fitness) * worldToScreenMult * -1, 250.0);
    scale(3.0);
  }
  
  // Check the 0% creature exists
  if (creature0th != null)
  {
    // Show the 0% creature in top left
    translate(0, 333);
    scale(0.333333);
    fill(140);
    rect(0, 0, 3000, 1000);
    fill(0);
    textSize(60);
    text(creatureCount + ". [[" + creature0th.DNA + "]]   " + creature0th.fitness + " metres", 950, 80);
    //translate((viewWidth / 2 - creature0th.fitness) * worldToScreenMult, -250.0);
    translate(600, 0.0);
    creature0th.Draw();
    translate(-600, 0.0);
    //translate((viewWidth / 2 - creature0th.fitness) * worldToScreenMult * -1, 250.0);
    scale(3.0);
  }
}

/////////////////////////////////////////////
//                                         //
//              MAIN PROGRAM               //
//                                         //
/////////////////////////////////////////////

// All the creatures to simulate
SquishyCreature[] creatures;
float[] distances;
int winningCreature = -1;

// Simulates a creature. Returns the distance that
// the creature travelled
float SimulateCreature(SquishyCreature creature, int numSteps, float timestep)
{
  // Loop over steps
  for (int i = 0; i < numSteps; i++)
    creature.PhysicsStep(timestep);
  
  // Fetch the list of nodes
  BouncyNode[] nodes = creature.nodes;
  
  // Add up the average position
  float averagePos = 0;
  int sumCount = 0;
  
  // Loop over each node
  for (int i = 0; i < nodes.length; i++)
  {
    // Check node exists
    if (nodes[i] != null)
    {
      // Don't add to average if NaN position
      if (!Float.isNaN(nodes[i].position.x) && !Float.isNaN(nodes[i].position.y) && nodes[i].position.x < 1000.0 && nodes[i].position.y < 1000.0)
      {
        averagePos += nodes[i].position.x;
        sumCount++;
      }
    }
  }
  
  // Average out the node position
  averagePos /= float(sumCount);
  
  // If sum count is zero, give zero score
  if (sumCount < 1)
    averagePos = 0;
  
  // Return the magnitude of the distance travelled
  return averagePos;
}

// Simulates a generation
void SimulateGeneration()
{
  // Reset pre-sim timer
  presimTimer = 0;
  
  // Batch is now in progress
  batchInProgress = true;
  
  // Time before
  float timeBefore = millis();
  
  // Simulate all creatures
  for (int i = 0; i < creatureCount; i++)
  {
    // Simulate creatures
    distances[i] = SimulateCreature(creatures[i], int(simTime / simTimestep), simTimestep);
    
    // Assign fitness
    creatures[i].fitness = distances[i];
  }
  
  // Find the highest distance
  float maxDistance = Float.NEGATIVE_INFINITY;
  float avgDistance = 0.0;
  SquishyCreature bestCreature = CreatureFromDNA(GenerateNewDNA(startGridWidth, startGridHeight));
  for (int i = 0; i < creatureCount; i++)
  {
    avgDistance += distances[i];
    if (distances[i] > maxDistance)
    {
      maxDistance = distances[i];
      bestCreature = creatures[i];
    }
  }
  avgDistance /= float(creatureCount);
  
  // Sort the list using Java's built-in array sort
  // Relies on SquishyCreature implementing Comparable<T>
  Arrays.sort(creatures);
  
  // Set significant percentile creatures
  //creatureA = bestCreature;
  creature0th = creatures[0];
  creature50th = creatures[creatureCount / 2 - 1];
  creature100th = creatures[creatureCount - 1];
  
  
  // Reset position
  for (int y = 0; y < startGridHeight + 1; y++)
    for (int x = 0; x < startGridWidth + 1; x++)
      if (creature0th.nodes[x + y * (startGridWidth + 1)] != null)
        creature0th.nodes[x + y * (startGridWidth + 1)].position = PVector.add(startPosition, new PVector(x * startCellSize, y * startCellSize));
  
  // Reset position
  for (int y = 0; y < startGridHeight + 1; y++)
    for (int x = 0; x < startGridWidth + 1; x++)
      if (creature50th.nodes[x + y * (startGridWidth + 1)] != null)
        creature50th.nodes[x + y * (startGridWidth + 1)].position = PVector.add(startPosition, new PVector(x * startCellSize, y * startCellSize));
  
  // Reset position
  for (int y = 0; y < startGridHeight + 1; y++)
    for (int x = 0; x < startGridWidth + 1; x++)
      if (creature100th.nodes[x + y * (startGridWidth + 1)] != null)
        creature100th.nodes[x + y * (startGridWidth + 1)].position = PVector.add(startPosition, new PVector(x * startCellSize, y * startCellSize));
  
  
  // Create a file to save
  String[] saveData = new String[creatureCount];
  for (int i = 0; i < saveData.length; i++)
    saveData[i] = (creatureCount - i) + ". [[" + creatures[i].DNA + "]] " + creatures[i].fitness + "\n";
  
  // Save the generation data
  if (saveGenerations)
    saveStrings("generation " + (generationCount + 1) + ".txt", saveData);
  
  
  // Kill 500 creatures. Using a gaussian/normal distribution, we will make it so
  // it will generally be the weakest creatures that are killed, but stronger creatures
  // still have a chance to die
  for (int i = 0; i < creatureCount / 2; i++)
  {
    // Get the number and scale it to our data set
    float pos = randomGaussian() * creatureCount;
    
    // If it's below zero, multiply by -1 to make it positive
    if (pos < 0) pos *= -1;
    
    // Make sure we kill a creature
    // I have no idea if this statement will work??
    while ((int(pos) >= creatureCount ? true : creatures[int(pos)] == null) || int(pos) >= creatureCount)
    {
      pos = randomGaussian() * creatureCount;
      if (pos < 0) pos *= -1;
    }
    
    // Remove this creature
    creatures[int(pos)] = null;
  }
  
  
  // Create a new generation
  ArrayList<SquishyCreature> newGeneration = new ArrayList<SquishyCreature>();
  
  // Loop over each surviving creature, and make them reproduce twice
  for (int i = 0; i < creatures.length; i++)
  {
    if (creatures[i] != null)
    {
      newGeneration.add(creatures[i].Copy());
      newGeneration.add(creatures[i].Copy());
    }
  }
  
  // Set the creatures list again
  //creatures = (String[])newGeneration.toArray();
  for (int i = 0; i < creatureCount; i++)
    creatures[i] = newGeneration.get(i);
  
  // Time after simulating
  float timeAfter = millis();
  
  // Announce time taken
  print("Generation " + generationCount + " took " + (timeAfter - timeBefore) + " ms. Best: " + bestCreature.fitness + " [[" + bestCreature.DNA + "]], average: " + avgDistance + "\n");
  
  // Increment gen count
  generationCount++;
  
  // Batch is no longer in process
  batchInProgress = false;
}

// Called on start
void setup()
{
  // Set window size
  size(1000, 1000);
  
  // Seed the number generator
  randomSeed(seed);
  
  // Initialise arrays
  creatures = new SquishyCreature[creatureCount];
  distances = new float[creatureCount];
  
  // Create the initial generation
  for (int i = 0; i < creatureCount; i++)
    creatures[i] = CreatureFromDNA(GenerateNewDNA(startGridWidth, startGridHeight));
  
  // If simulating automatically, begin
  if (maxAutoGens > 0)
    simulating = true;
  
  // Setup test creatures
  //creature0th = CreatureFromDNA("_____#_______________________#_____#111111111111111111111111111110111111");
  creature0th = CreatureFromDNA(GenerateNewDNA(startGridWidth, startGridHeight));
  creature50th = CreatureFromDNA(GenerateNewDNA(startGridWidth, startGridHeight));
  creature100th = CreatureFromDNA(GenerateNewDNA(startGridWidth, startGridHeight));
  //creature0th = CreatureFromDNA("#########111111111");
}

// Called on keypress
void keyPressed()
{
  // Change view type
  //creatureViewType++;
  
  // Simulate another generation
  //SimulateGeneration();
  
  // If simulating, stop
  if (simulating) simulating = false;
  else simulating = true;
  
  // Say whether we're simulating
  print(simulating == true ? "Simulation started.\n" : "Simulation stopped.\n");
}

// Called every frame
void draw()
{
  // Stop if we've exceeded max number of auto-simulated generations
  if (generationCount >= maxAutoGens && maxAutoGens > 0)
  {
    simulating = false;
    exit();
  }
  
  // Increment presim timer
  presimTimer += 1.0 / 60.0;
  
  // Set a background
  background(85, 185, 225);
  
  // Keep creature view type in bounds
  if (creatureViewType > 2)
    creatureViewType = 0;
  
  // If simulating, do simulations
  if (simulating)
    SimulateGeneration();
  
  // Draw the menu
  drawUI();
  
  // Do physics when not simulating
  if (!simulating)
  {
    // Time step. Making this lower (1/300 is the base used, 1/450 is nice, 1/600 is even better) will
    // improve smoothness but reduce simulation speed
    float timestep = 1.0 / 600.0;
    
    // How much to speed up time. 10 makes it 10x faster. Note: with timestep set to 1/600,
    // 10x speed is realtime 60 fps
    int timewarp = 10;
    
    // Warp time
    for (int i = 0; i < timewarp; i++)
    {
      // Simulate each creature
      if (creature0th != null) creature0th.PhysicsStep(timestep);
      if (creature50th != null) creature50th.PhysicsStep(timestep);
      if (creature100th != null) creature100th.PhysicsStep(timestep);
    }
  }
}
