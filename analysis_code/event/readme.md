# CLAS12 Event Class
This class was designed to be a useful data structure to hold the CLAS12 detector information.  Each detector is represented in a map where the index of the particle in `REC::Particle` defines the map key.

### Simple Use
You can convert an entire `HipoDataEvent` using the `EventConverter.convert` method.

```groovy
def reader = HipoDataSource()
reader.open(filename)

while(reader.hasEvent()){
    def dataEvent = reader.getNextEvent()

    // This is where the conversion is done.
    def event = EventConverter.convert(dataEvent)
}

```

### Advanced Use
If you only need to check banks other than `REC::Particle` in a small fraction of your events, you can call the bank converters directly to increase the execution speed and reduce the memory usage of each event.  

```groovy
def reader = HipoDataSource()
reader.open(filename)

while(reader.hasEvent()){
    def dataEvent = reader.getNextEvent()

    def event = new Event()
    EventConverter.convertScalar(dataEvent, event)
    EventConverter.convertPart(dataEvent, event)

    // Check something in REC::Particle and call other banks if needed. 

}

```

### Event Structure
The `Event` class contains detector information indexed by the index of the particle in `REC::Particle`.  Let's use the number of photoelectrons in the Cherenkov counter as an example.  Here, we select event builder electrons and get the number of photoelectrons.  You should remember to check that the particle exists in the detector before trying to access it (or return a null as default and handle that later).  In order to check, a HashSet has been added for each detector that contains the indices of all particles that left signal in that detector.  In this example, `event.cherenkov_status` is used. 

```groovy

(0 ..< event.npart).findAll{ index ->
   event.pid.get(index) == 11 
}.collect{ index ->
   if (event.cherenkov_status.contains(index)) { event.nphe.get(index) }
}

```

Please note that in groovy you can use the explicit map get as above, or the shorter notation `event.pid[index]`.  I have opted to use the first notation simply to reinforce that the underlying data structure is a map.

### Drift Chamber Structure
The drift chamber structure is slightly more complicated.  I have introduced the class `DCHit` to keep track of all the hits for each particle.  To take a look at all the hits left by a particle with `REC::Particle` index you can do:

```groovy

if (event.dc1_status.contains(index)){
   def hits = event.dc1.get(index)

   // the number of hits
   def nhits = hits.size()

   // all of the x-y coordinates 
   def xcoords = hits*.x
   def ycoords = hits*.y

   // x-y coordinates for hits in
   // the last 6 layers
   def coords = hits.findAll{ hit -> hit.layer > 6 }.collect{ hit -> [hit.x, hit.y]}

   // the entire first hit
   def hit = hits.getAt(0)
   def hit = hits[0]

   // the hit in layer 12, if it exists
   def hit = hits.find{it.layer==12}
}

```



