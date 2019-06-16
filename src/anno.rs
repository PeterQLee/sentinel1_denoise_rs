struct SwathBounds {
    firstAzimuthLine:usize,
    lastAzimuthLine:usize,
    firstRangeSample:usize,
    lastRangeSample:usize
}

impl SwathBounds {
    //fn merge_parse_list(reader:&mut Reader<&[u8]>)
    fn swath_parse_func(reader:&mut Reader<&[u8]>) -> SwathBounds {
        //swathmerging/swathmergelist/swathboundslist/(swathbounds)

        let mut bounds = SwathBounds{
            firstAzimuthLine:0,
            lastAzimuthLine:0,
            firstRangeSample:0,
            lastRangeSample:0
        };
        enum ExtractState {FirstAz, LastAz, FirstRg, LastRg, None};
        loop {
            match reader.read_event(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"swathBounds" => {state = ExtractState::None},
                        b"firstAzimuthLine" => {state = ExtractState::FirstAz},
                        b"lastAzimuthLine" => {state = ExtractState::LastAz},
                        b"firstRangeSample" => {state = ExtractState::FirstRg},
                        b"lastRangeSample" => {state = ExtractState::LastRg}
                        _ => {},
                    }
                }
                Ok(Event::Text(ref e)) =>{
                    match state {
                        match ExtractState::FirstAz => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => bounds.firstAzimuthLine = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                            
                        }
                        match ExtractState::LastAz => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => bounds.lastAzimuthLine = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        }
                        match ExtractState::FirstRg => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => bounds.firstRangeLine = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                            
                        }
                        match ExtractState::LastRg => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => bounds.lastRangeLine = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"swathBounds"=>{break;},
                        b"firstAzimuthLine" => {state = ExtractState::None},
                        b"lastAzimuthLine" => {state = ExtractState::None},
                        b"firstRangeSample" => {state = ExtractState::None},
                        b"lastRangeSample" => {state = ExtractState::None},
                        _ => {}
                    }
                }
                _ => {},
            }
        }
        return bounds;
    }
    
}


struct Parameters {
    numswath:usize
    numAzimuth:Vec<usize>,
    numBurst:Vec<usize>,
    swathBoundList:Vec<SwathBounds>,
}
