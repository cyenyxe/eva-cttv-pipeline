/*
 *  Copyright 2018 EMBL - European Bioinformatics Institute
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
package uk.ac.ebi.eva.clinvar;

import org.junit.Test;

import java.io.InputStream;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;

public class XmlClinVarReaderTest {

    private static final int CLINVAR_RECORDS_IN_TEST_FILE = 10;

    @Test
    public void testCall() throws Exception {
        InputStream clinvarExampleInputStream = this.getClass().getResourceAsStream("/ClinvarExample.xml");
        // the capacity of the queue is the number of records in the test file plus one, because the reader
        // adds an object at the end to indicate it has finished reading
        ArrayBlockingQueue<String> queue = new ArrayBlockingQueue<>(CLINVAR_RECORDS_IN_TEST_FILE + 1);
        XmlClinVarReader reader = new XmlClinVarReader(clinvarExampleInputStream, queue);

        int clinvarRecordsRead = reader.call();

        // collect the queue objects into a list because it is more convenient to assert the results
        List<String> clinvarRecords = queue.stream().collect(Collectors.toList());

        assertEquals(CLINVAR_RECORDS_IN_TEST_FILE, clinvarRecordsRead);
        assertEquals(CLINVAR_RECORDS_IN_TEST_FILE + 1, clinvarRecords.size());

        // the last record in the queue is the one that the reader adds to indicate it has finished
        assertEquals(XmlClinVarReader.FINISHED, clinvarRecords.get(CLINVAR_RECORDS_IN_TEST_FILE));
    }
}